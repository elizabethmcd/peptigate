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

# Get all genome assembly files ending in "*.fa" from input directory
genome_files = {f.stem: str(f) for f in INPUT_DIR.glob("*.fa")}
print(genome_files)
genome_name = list(genome_files.keys())
print(genome_name)

################################################################################
## Prodigal protein prediction from MAGs
################################################################################

# All combined ORFs AAs from input MAGs
ORFS_AMINO_ACIDS = OUTPUT_DIR / "prodigal" / "all_mag_proteins.faa"

rule prodigal_predictions:
    """
    Run Prodigal to predict protein-coding genes and translate them into proteins.
    """
    input:
        fasta=lambda wildcards: genome_files[wildcards.genome_name]  # Input genome FASTA from samplesheet
    output:
        gbk=OUTPUT_DIR / "prodigal" / "{genome_name}.gbk",
        faa=OUTPUT_DIR / "prodigal" / "{genome_name}.faa",
        fna=OUTPUT_DIR / "prodigal" / "{genome_name}.fna",
    conda:
        "envs/prodigal.yml"
    shell:
        """
        prodigal -i {input.fasta} -o {output.gbk} -a {output.faa} -d {output.fna}
        """

rule combine_prodigal_outputs:
    """
    Combine all Prodigal faa outputs into a single file for downstream use.
    """
    input:
        expand(OUTPUT_DIR / "prodigal" / "{genome_name}.faa", genome_name=genome_name)
    output:
        ORFS_AMINO_ACIDS
    shell:
        """
        cat {input} > {output}
        """

################################################################################
## smORF prediction
################################################################################


rule smorfinder_prediction:
    """
    Use smorfinder to identify small ORFs (smORFs) in the input protein sequences.
    """
    input:
        fasta=lambda wildcards: genome_files[wildcards.genome_name]
    output:
        gff=OUTPUT_DIR / "smORF" / "smorfinder_predictions.gff",
        peptide_faa=OUTPUT_DIR / "smORF" / "smorf_peptides.faa",
        ffn=OUTPUT_DIR / "smORF" / "smorf_nucleotide_sequences.ffn",
        tsv=OUTPUT_DIR / "smORF" / "smorf_summary.tsv",
    conda:
        "envs/smorfinder.yml"
    shell:
        """
        smorf single {input.fasta} -o {wildcards.genome_name}_smorf_output
        mv {wildcards.genome_name}_smorf_output.gff {output.gff}
        mv {wildcards.genome_name}_smorf_output.faa {output.peptide_faa}
        mv {wildcards.genome_name}_smorf_output.ffn {output.ffn}
        mv {wildcards.genome_name}_smorf_output.tsv {output.tsv}
        """


rule convert_smorf_peptide_faa_to_tsv:
    input:
        peptide_faa=rules.smorfinder_prediction.output.peptide_faa,
    output:
        tsv=OUTPUT_DIR / "smORF" / "smorf_peptides.tsv",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit fx2tab --only-id {input} -o {output}
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
        parent_faa=OUTPUT_DIR / "cleavage" / "nlpprecursor" / "nlpprecursor_peptide_parents.faa",
        peptide_faa=OUTPUT_DIR / "cleavage" / "nlpprecursor" / "nlpprecursor_peptides.faa",
        tsv=OUTPUT_DIR / "cleavage" / "nlpprecursor" / "nlpprecursor_predictions.tsv",
    params:
        modelsdir=INPUT_DIR / "models" / "nlpprecursor" / "models/",
    conda:
        "envs/nlpprecursor.yml"
    shell:
        """
        python scripts/run_nlpprecursor.py \
            --models_dir {params.modelsdir} \
            --protein_fasta_file {input.faa} \
            --proteins_output_file {output.parent_faa} \
            --protein_peptides_output_file {output.peptide_faa} \
            --predictions_output_file {output.tsv}
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
        json=rules.deeppeptide.output.json,
        faa=rules.remove_stop_codon_asterisk_from_transdecoder_ORFs.output.faa,
    output:
        parent_faa=OUTPUT_DIR / "cleavage" / "deeppeptide" / "deeppeptide_peptide_parents.faa",
        peptide_faa=OUTPUT_DIR / "cleavage" / "deeppeptide" / "deeppeptide_peptides.faa",
        tsv=OUTPUT_DIR / "cleavage" / "deeppeptide" / "deeppeptide_predictions.tsv",
    conda:
        "envs/biopython.yml"
    shell:
        """
        python scripts/extract_deeppeptide_sequences.py \
            --json_file {input.json} \
            --protein_fasta_file {input.faa} \
            --proteins_output_file {output.parent_faa} \
            --protein_peptides_output_file {output.peptide_faa} \
            --predictions_output_file {output.tsv}
        """


################################################################################
## Combine peptide predictions
################################################################################


rule combine_peptide_faa_predictions:
    input:
        nlpprecursor=rules.nlpprecursor.output.peptide_faa,
        deeppeptide=rules.extract_deeppeptide_sequences.output.peptide_faa,
        smorfinder=rules.smorfinder_prediction.output.peptide_faa  # Add smorfinder output
    output:
        peptide_faa=OUTPUT_DIR / "predictions" / "peptides.faa",
    shell:
        """
        cat {input} > {output.peptide_faa}
        """


rule convert_peptide_faa_to_tsv:
    input:
        peptide_faa=rules.combine_peptide_faa_predictions.output.peptide_faa,
    output:
        tsv=OUTPUT_DIR / "predictions" / "peptides_faa.tsv",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit fx2tab --only-id {input} -o {output}
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
        db=INPUT_DIR / "databases" / "peptipedia.dmnd",
    params:
        dbprefix=INPUT_DIR / "databases" / "peptipedia",
    conda:
        "envs/diamond.yml"
    shell:
        """
        diamond makedb --in {input.db} -d {params.dbprefix}
        """


rule diamond_blastp_peptide_predictions_against_peptipedia_database:
    input:
        db=rules.make_diamond_db_from_peptipedia_database.output.db,
        peptide_faa=rules.combine_peptide_faa_predictions.output.peptide_faa,
    output:
        tsv=OUTPUT_DIR / "annotation" / "peptipedia" / "blastp_matches.tsv",
    params:
        dbprefix=INPUT_DIR / "databases" / "peptipedia",
    conda:
        "envs/diamond.yml"
    shell:
        """
        diamond blastp -d {params.dbprefix} -q {input.peptide_faa} -o {output.tsv} --header simple \
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
        peptide_faa=rules.combine_peptide_faa_predictions.output.peptide_faa,
    output:
        tsv=OUTPUT_DIR / "annotation" / "deepsig" / "deepsig_predictions.tsv",
    conda:
        "envs/deepsig.yml"
    shell:
        """
        deepsig -f {input} -o {output}.tmp -k euk
        python scripts/add_header_to_deepsig_tsv.py {output}.tmp {output}
        """


rule characterize_peptides:
    input:
        peptide_faa=rules.combine_peptide_faa_predictions.output.peptide_faa,
    output:
        tsv=OUTPUT_DIR / "annotation" / "characteristics" / "peptide_characteristics.tsv",
    conda:
        "envs/peptides.yml"
    shell:
        """
        python scripts/characterize_peptides.py {input.peptide_faa} {output.tsv}
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
        peptide_faa=rules.combine_peptide_faa_predictions.output.peptide_faa,
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
            --input_fasta {input.peptide_faa} \
            --model_folder {params.modelsdir}/{wildcards.autopeptideml_model_name}_1/ensemble \
            --model_name {wildcards.autopeptideml_model_name} \
            --output_tsv {output.tsv}
        """


## TER TODO: remove ffn and adapt script; remove plmutils as well
rule combine_peptide_predictions:
    input:
        nlpprecursor=rules.nlpprecursor.output.tsv,
        deeppeptide=rules.extract_deeppeptide_sequences.output.tsv,
        sorf=rules.convert_smorf_peptide_faa_to_tsv.output.tsv,
        faa_tab=rules.convert_peptide_faa_to_tsv.output.tsv,
    output:
        tsv=OUTPUT_DIR / "predictions" / "peptide_predictions.tsv",
    conda:
        "envs/tidyverse.yml"
    shell:
        """
        Rscript scripts/combine_peptide_predictions_protein_input.R \
            --nlpprecursor_path {input.nlpprecursor} \
            --deeppeptide_path {input.deeppeptide} \
            --sorf_path {input.sorf} \
            --faa_tab_path {input.faa_tab} \
            --output_predictions_path {output.tsv}
        """


rule combine_peptide_annotations:
    input:
        autopeptideml=expand(
            rules.run_autopeptideml.output.tsv, autopeptideml_model_name=AUTOPEPTIDEML_MODEL_NAMES
        ),
        deepsig=rules.run_deepsig.output.tsv,
        peptipedia=rules.diamond_blastp_peptide_predictions_against_peptipedia_database.output.tsv,
        characteristics=rules.characterize_peptides.output.tsv,
    output:
        tsv=OUTPUT_DIR / "predictions" / "peptide_annotations.tsv",
    params:
        autopeptidemldir=OUTPUT_DIR / "annotation" / "autopeptideml/",
    conda:
        "envs/tidyverse.yml"
    shell:
        """
        Rscript scripts/combine_peptide_annotations.R \
            --autopeptideml_dir {params.autopeptidemldir} \
            --deepsig_path {input.deepsig} \
            --peptipedia_path {input.peptipedia} \
            --characteristics_path {input.characteristics} \
            --output_annotations_path {output.tsv}
        """


################################################################################
## Target rule all
################################################################################


rule all:
    default_target: True
    input:
        rules.combine_peptide_predictions.output.tsv,
        rules.combine_peptide_annotations.output.tsv,
