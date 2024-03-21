import os
from pathlib import Path

################################################################################
## Configuration
################################################################################

# Retrieves the absolute path of the directory snakemake is launched in.
# Used by DeepPeptide to simplify output file paths.
WORKING_DIRPATH = Path(os.getcwd())


# Default pipeline configuration parameters are in the config file.
# If you create a new yml file and use the --configfile flag,
# options in that new file overwrite the defaults.
configfile: "./config.yml"


INPUT_DIR = Path(config["input_dir"])
OUTPUT_DIR = Path(config["output_dir"])

ORFS_AMINO_ACIDS = Path(config["orfs_amino_acids"])
ORFS_NUCLEOTIDES = Path(config["orfs_nucleotides"])
CONTIGS_SHORTER = Path(config["contigs_shorter_than_r2t_minimum_length"])
CONTIGS_LONGER = Path(config["contigs_longer_than_r2t_minimum_length"])
PLMUTILS_MODEL_DIR = Path(config["plmutils_model_dir"])

################################################################################
## sORF prediction
################################################################################


rule combine_contigs:
    """
    By default we assume that files provided to this pipeline are reads2transcriptome outputs.
    Reads2transcriptome outputs two files that contain contigs.
    The first we refer to as contigs_shorter_than_r2t_minimum (or CONTIGS_SHORTER),
    which are transcripts output by assemblers that did not meet the r2t runs minimum contig length.
    The second we refer to as contigs_longer_than_r2t_minimum (or CONTIGS_LONGER) and are assembled
    transcripts that passed the isoform clustering and decontamination steps of r2t.
    We no longer need these pools of transcripts differentiated, so we combine them in this rule. 
    If your input transcriptome only has one file of assembled contigs, sequences only need to be
    supplied in contigs_longer_than_r2t_minimum_length.
    contigs_shorter_than_r2t_minimum_length can be an empty file.
    """
    input:
        contigs_shorter=CONTIGS_SHORTER,
        contigs_longer=CONTIGS_LONGER,
    output:
        all_contigs=OUTPUT_DIR / "sORF" / "contigs" / "all_input_contigs.fa",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        cat {input.contigs_shorter} {input.contigs_longer} > {output.all_contigs}
        """


rule get_coding_contig_names:
    """
    Extract amino acid contig names and remove everything after the first period, 
    which are isoform labels.
    This file will be used to select all contigs that DID NOT have a transdecoder-detected ORFs. 
    """
    input:
        ORFS_AMINO_ACIDS,
    output:
        names=OUTPUT_DIR / "sORF" / "contigs" / "orfs_amino_acid_names.txt",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit seq -n {input} | sed 's/[.].*$//' > {output}
        """


rule filter_contigs_to_no_predicted_ORF:
    """
    The r2t pipeline runs transdecoder to predict open reading frames (ORFs) from transcripts.
    By default, only ORFs that are longer than 100 amino acids are kept by transdecoder.
    The peptigate pipeline predicts peptides that are 100 amino acids or shorter.
    This rule eliminates transcripts that contained a transdecoder-predicted ORF.
    It keeps all other transcripts, regardless of length, to investigate the presence of an sORF
    later in the pipeline.
    """
    input:
        fa=rules.combine_contigs.output.all_contigs,
        names=rules.get_coding_contig_names.output.names,
    output:
        fa=OUTPUT_DIR / "sORF" / "contigs" / "contigs_with_no_transdecoder_predicted_orf.fa",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit grep -v -f {input.names} {input.fa} -o {output.fa}
        """


rule plmutils_translate:
    """
    This rule takes input nucleotide transcripts, detects the longest open reading frame, and
    translates it into amino acid sequences.
    """
    input:
        rules.filter_contigs_to_no_predicted_ORF.output.fa,
    output:
        faa=OUTPUT_DIR / "sORF" / "plmutils" / "translated_contigs.faa",
    conda:
        "envs/plmutils.yml"
    shell:
        """
        plmutils translate --longest-only --output-filepath {output} {input}
        """


rule length_filter_plmutils_translate_output:
    input:
        rules.plmutils_translate.output.faa,
    output:
        faa=OUTPUT_DIR / "sORF" / "plmutils" / "translated_contigs_filtered.faa",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit seq --max-len 100 -o {output} {input}
        """


rule plmutils_embed:
    """
    This rule embeds amino acid sequences produced by plmutils translated into the embedding space
    of a protein large language model.
    For now, plmutils only supports ESM.
    We use the 8M model because it is fast and there is some evidence that it works well for other
    peptide prediction tasks (https://www.biorxiv.org/content/10.1101/2023.11.13.566825v2).
    The parameter --layer-ind -1 means to extract the embedding from the last layer of the model.
    """
    input:
        rules.length_filter_plmutils_translate_output.output.faa,
    output:
        npy=OUTPUT_DIR / "sORF" / "plmutils" / "embedded_contigs_filtered.npy",
    conda:
        "envs/plmutils.yml"
    shell:
        """
        plmutils embed --model-name esm2_t6_8M_UR50D \
            --layer-ind -1 \
            --output-filepath {output.npy} \
            {input}
        """


rule plmutils_predict:
    input:
        embeddings=rules.plmutils_embed.output.npy,
        faa=rules.length_filter_plmutils_translate_output.output.faa,
        model=PLMUTILS_MODEL_DIR,
    output:
        csv=OUTPUT_DIR / "sORF" / "plmutils" / "predictions.csv",
    conda:
        "envs/plmutils.yml"
    shell:
        """
        plmutils predict --model-dirpath {input.model} \
            --embeddings-filepath {input.embeddings} \
            --fasta-filepath {input.faa} \
            --output-filepath {output.csv}
        """


rule extract_plmutils_predicted_peptides:
    input:
        csv=rules.plmutils_predict.output.csv,
        faa=rules.length_filter_plmutils_translate_output.output.faa,
    output:
        names=OUTPUT_DIR / "sORF" / "plmutils" / "peptide_names.faa",
        faa=OUTPUT_DIR / "sORF" / "plmutils" / "peptides.faa",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        grep "positive" {input.csv} | cut -d, -f1 > {output.names} 
        seqkit grep -f {output.names} {input.faa} -o {output.faa}
        """


################################################################################
## cleavage prediction
################################################################################


rule remove_stop_codon_asterisk_from_transdecoder_ORFs:
    input:
        ORFS_AMINO_ACIDS,
    output:
        faa=OUTPUT_DIR / "cleavage" / "preprocessing" / "noasterisk.faa",
    shell:
        """
        sed '/^[^>]/s/\*//g' {input} > {output}
        """


# Ribosomally synthesized and post-translationally modified peptide prediction


rule download_nlpprecursor_models:
    output:
        tar=INPUT_DIR / "models" / "nlpprecursor" / "nlpprecursor_models.tar.gz",
        model=INPUT_DIR / "models" / "nlpprecursor" / "models" / "annotation" / "model.p",
    params:
        outdir=INPUT_DIR / "models" / "nlpprecursor",
    shell:
        """
        curl -JLo {output.tar} https://github.com/magarveylab/NLPPrecursor/releases/download/1.0/nlpprecursor_models.tar.gz
        tar xf {output.tar} -C {params.outdir} 
        """


rule filter_protein_sequences_with_nonstandard_amino_acids:
    """
    The NLPPrecursor tool only uses the 20 standard amino acids.
    This script removes sequences with nonstandard amino acids in the sequence.
    """
    input:
        faa=rules.remove_stop_codon_asterisk_from_transdecoder_ORFs.output.faa,
    output:
        faa=OUTPUT_DIR / "cleavage" / "preprocessing" / "noasterisk_nononstandardaa.faa",
    conda:
        "envs/nlpprecursor.yml"
    shell:
        """
        python scripts/filter_protein_sequences_with_nonstandard_amino_acids.py --input {input} --output {output}
        """


rule nlpprecursor:
    """
    The nlpprecursor tool is part of [DeepRiPP](https://doi.org/10.1073/pnas.1901493116)
    that predicts ribosomally synthesized postranslationally modified peptides, 
    a subclass of cleavage peptides.
    Unlike many tools in the RiPP prediction space, "NLPPrecursor identifies RiPPs independent of 
    genomic context and neighboring biosynthetic genes."
    Note the paper suggests, "the precursor cleavage algorithm predicted N-terminal cleavage sites 
    with 90% accuracy, when considering cleavage points ±5 amino acids from the true prediction 
    site, a range within which all possible complete chemical structures can be elaborated in 
    silico by combinatorial structure prediction."    
    From the DeepRiPP paper supplement: 
    "Protein sequences of open reading frames are used as input. The output of the model consists 
     of a classification of each ORF as either a precursor peptide (further subclassified according 
    to RiPP family), or a non-precursor peptide. A total of 14 classes are identified (n_class)."
    """
    input:
        faa=rules.filter_protein_sequences_with_nonstandard_amino_acids.output.faa,
        model=rules.download_nlpprecursor_models.output.model,
    output:
        tsv=OUTPUT_DIR / "cleavage" / "nlpprecursor" / "nlpprecursor_predictions.tsv",
        peptide=OUTPUT_DIR / "cleavage" / "nlpprecursor" / "nlpprecursor_peptides.fasta",
    params:
        modelsdir=INPUT_DIR / "models" / "nlpprecursor" / "models/",
    conda:
        "envs/nlpprecursor.yml"
    shell:
        """
        python scripts/run_nlpprecursor.py {params.modelsdir} {input.faa} {output.tsv} {output.peptide}
        """


# General Cleavage peptide prediction


rule clone_deeppeptide:
    output:
        src=touch("cloned_repositories/DeepPeptide/deeppeptide_cloned.txt"),
    shell:
        """
        # only clone the repo if it isn't already present
        if [ ! -d cloned_repositories/DeepPeptide ]; then
            git clone https://github.com/fteufel/DeepPeptide.git cloned_repositories/DeepPeptide
        fi
        cd cloned_repositories/DeepPeptide
        git checkout 2657f5dca38e6417c65da5913c1974ed932746e3
        """


rule deeppeptide:
    input:
        src=rules.clone_deeppeptide.output.src,
        faa=rules.remove_stop_codon_asterisk_from_transdecoder_ORFs.output.faa,
    output:
        json=OUTPUT_DIR / "cleavage" / "deeppeptide" / "peptide_predictions.json",
    conda:
        "envs/deeppeptide.yml"
    params:
        outdir=OUTPUT_DIR / "cleavage" / "deeppeptide",
    shell:
        """
        cd cloned_repositories/DeepPeptide/predictor && python3 predict.py --fastafile {WORKING_DIRPATH}/{input.faa} --output_dir {params.outdir} --output_fmt json
        mv {params.outdir}/* {WORKING_DIRPATH}/{params.outdir}/
        """


rule extract_deeppeptide_sequences:
    """
    DeepPeptide outputs a json file of peptide predictions and locations,
    but does not output the sequences themselves.
    This step parses the JSON file and the protein FASTA from which peptides were predicted.
    It outputs the propeptide (full ORF, uncleaved) and the predicted peptide sequence (cleaved)
    in FASTA format.
    """
    input:
        faa=rules.remove_stop_codon_asterisk_from_transdecoder_ORFs.output.faa,
        json=rules.deeppeptide.output.json,
    output:
        propeptide=OUTPUT_DIR / "cleavage" / "deeppeptide" / "propeptides.faa",
        peptide=OUTPUT_DIR / "cleavage" / "deeppeptide" / "peptides.faa",
        tsv=OUTPUT_DIR / "cleavage" / "deeppeptide" / "predictions.tsv",
    conda:
        "envs/biopython.yml"
    shell:
        """
        python scripts/extract_deeppeptide_sequences.py {input.json} {input.faa} {output.propeptide} {output.peptide} {output.tsv}
        """


################################################################################
## Combine peptide predictions
################################################################################


rule combine_peptide_predictions:
    input:
        sorf=rules.extract_plmutils_predicted_peptides.output.faa,
        nlpprecursor=rules.nlpprecursor.output.peptide,
        deeppeptide=rules.extract_deeppeptide_sequences.output.peptide,
    output:
        peptide=OUTPUT_DIR / "annotation" / "combined_peptide_predictions" / "peptides.faa",
    shell:
        """
        cat {input} > {output.peptide}
        """


################################################################################
## Compare against known peptides
################################################################################


rule download_peptipedia_database:
    """
    The peptipedia database includes sequences from 66 peptide databases.
    It was last updated in 01/2024.
    We selected this database because it has collected the most peptides in a single location and 
    done some quality control on those sequences.
    """
    output:
        db=INPUT_DIR / "databases" / "peptipedia.fasta.gz",
    shell:
        """
        curl -JLo {output} https://osf.io/dzycu/download 
        """


rule make_diamond_db_from_peptipedia_database:
    input:
        db=rules.download_peptipedia_database.output.db,
    output:
        db=OUTPUT_DIR / "annotation" / "peptipedia" / "0_diamond_db" / "peptipedia.dmnd",
    params:
        dbprefix=OUTPUT_DIR / "annotation" / "peptipedia" / "0_diamond_db" / "peptipedia",
    conda:
        "envs/diamond.yml"
    shell:
        """
        diamond makedb --in {input.db} -d {params.dbprefix}
        """


rule diamond_blastp_peptide_predictions_against_peptipedia_database:
    input:
        db=rules.make_diamond_db_from_peptipedia_database.output.db,
        peptide=rules.combine_peptide_predictions.output.peptide,
    output:
        tsv=OUTPUT_DIR / "annotation" / "peptipedia" / "1_blastp" / "matches.tsv",
    params:
        dbprefix=OUTPUT_DIR / "annotation" / "peptipedia" / "0_diamond_db" / "peptipedia",
    conda:
        "envs/diamond.yml"
    shell:
        """
        diamond blastp -d {params.dbprefix} -q {input.peptide} -o {output.tsv} --header simple \
         --outfmt 6 qseqid sseqid full_sseq pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore
        """


################################################################################
## Charaterize & annotate predicted peptide sequences
################################################################################


rule run_deepsig:
    """
    This rule uses deepsig to predict signal peptides in proteins using deep learning.
    """
    input:
        peptide=rules.combine_peptide_predictions.output.peptide,
    output:
        tsv=OUTPUT_DIR / "annotation" / "deepsig" / "deepsig.tsv",
    conda:
        "envs/deepsig.yml"
    shell:
        """
        deepsig -f {input} -o {output}.tmp -k euk
        python scripts/add_header_to_deepsig_tsv.py {output}.tmp {output}
        """


rule characterize_peptides:
    input:
        peptide=rules.combine_peptide_predictions.output.peptide,
    output:
        tsv=OUTPUT_DIR / "annotation" / "characteristics" / "peptide_characteristics.tsv",
    conda:
        "envs/peptides.yml"
    shell:
        """
        python scripts/characterize_peptides.py {input.peptide} {output.tsv}
        """


AUTOPEPTIDEML_MODEL_NAMES = [
    "AB",
    "ACE",
    "ACP",
    "AF",
    "AMAP",
    "AMP",
    "AOX",
    "APP",
    "AV",
    "BBP",
    "DPPIV",
    "MRSA",
    "Neuro",
    "QS",
    "TOX",
    "TTCA",
]


rule download_autopeptideml_models:
    output:
        # TER TODO: the authors of autopeptideml sent me these models.
        # They said they're working on uploading them.
        # Once they're available, I need to add a rule to download them and update the input here
        # to be the rules syntax
        # https://github.com/IBM/AutoPeptideML/issues/6#issuecomment-1989494559
        model=INPUT_DIR
        / "models"
        / "autopeptideml"
        / "HPO_NegSearch_HP"
        / "{autopeptideml_model_name}_1"
        / "apml_config.json",
    shell:
        """
        touch {output} # will become a curl command or something, depending on how models are packaged
        """


rule run_autopeptideml:
    """
    AutoPeptideML predicts the bioactivity of a peptide based on user-supplied models.
    The tool is a binary classifier, so each bioactivty has it's own model.
    As defined by AUTOPEPTIDEML_MODEL_NAMES, we use models trained in the autopeptideml preprint.
    The abbreviations are AB: Antibacterial; ACE: ACE inhibitor; ACP: Anticancer; AF: Antifungal;
    AMAP: Antimalarial; AMP: Antimicrobial; AOX: Antioxidant; APP: Antiparasitic; AV: Antiviral; 
    BBB: Brain-blood barrier crossing; DPPIV: DPPIV inhibitor; MRSA: Anti-MRSA; NP: Neuropeptide; 
    QS: Quorum sensing; TOX: Toxic; TTCA: Tumor T-cell antigens.
    
    The script below only predicts the bioactive classification against these models.
    However, autopeptideml was built to train new binary classifiers and peptipedia contains a lot
    of labelled peptides, so one could develop new models if the ones included above are 
    insufficient.
    """
    input:
        peptide=rules.combine_peptide_predictions.output.peptide,
        model=rules.download_autopeptideml_models.output.model,
    output:
        tsv=OUTPUT_DIR
        / "annotation"
        / "autopeptideml"
        / "autopeptideml_{autopeptideml_model_name}.tsv",
    params:
        modelsdir=INPUT_DIR / "models" / "autopeptideml" / "HPO_NegSearch_HP/",
    conda:
        "envs/autopeptideml.yml"
    shell:
        """
        python scripts/run_autopeptideml.py \
            --input_fasta {input.peptide} \
            --model_folder {params.modelsdir}/{wildcards.autopeptideml_model_name}_1/ensemble \
            --model_name {wildcards.autopeptideml_model_name} \
            --output_tsv {output.tsv}
        """


rule combine_peptide_annotations:
    input:
        nlpprecursor=rules.nlpprecursor.output.tsv,
        deeppeptide=rules.extract_deeppeptide_sequences.output.tsv,
        plmutils=rules.plmutils_predict.output.csv,
        autopeptideml=expand(
            rules.run_autopeptideml.output.tsv, autopeptideml_model_name=AUTOPEPTIDEML_MODEL_NAMES
        ),
        deepsig=rules.run_deepsig.output.tsv,
        peptipedia=rules.diamond_blastp_peptide_predictions_against_peptipedia_database.output.tsv,
        characteristics=rules.characterize_peptides.output.tsv,
    output:
        tsv1=OUTPUT_DIR / "annotation" / "peptide_predictions.tsv",
        tsv2=OUTPUT_DIR / "annotation" / "peptide_annotations.tsv",
    params:
        autopeptidemldir=OUTPUT_DIR / "annotation" / "autopeptideml/",
    conda:
        "envs/tidyverse.yml"
    shell:
        """
        Rscript scripts/combine_peptide_annotations.R \
            --nlpprecursor_path {input.nlpprecursor} \
            --deeppeptide_path {input.deeppeptide} \
            --plmutils_path {input.plmutils} \
            --autopeptideml_dir {params.autopeptidemldir} \
            --deepsig_path {input.deepsig} \
            --peptipedia_path {input.peptipedia} \
            --characteristics_path {input.characteristics} \
            --output_predictions_path {output.tsv1} \
            --output_annotations_path {output.tsv2}
        """


################################################################################
## Target rule all
################################################################################


rule all:
    default_target: True
    input:
        rules.combine_peptide_annotations.output.tsv1,


rule predict_sORF:
    """
    Defines a target rule for sORF prediction so a user can run only sORF prediction.
    snakemake predict_sORF --software-deployment-method conda -j 8 
    """
    input:
        sorf=rules.extract_plmutils_predicted_peptides.output.faa,


rule predict_cleavage:
    """
    Defines a target rule for cleavage prediction so a user can run only cleavage prediction.
    snakemake predict_cleavage --software-deployment-method conda -j 8 
    """
    input:
        rules.nlpprecursor.output.peptide,
        rules.extract_deeppeptide_sequences.output.peptide,
