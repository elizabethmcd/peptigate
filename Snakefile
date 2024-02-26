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

SHORT_CONTIGS = Path(config["short_contigs"])
ORFS_AMINO_ACIDS = Path(config["orfs_amino_acids"])
ORFS_NUCLEOTIDES = Path(config["orfs_nucleotides"])
ALL_CONTIGS = Path(config["all_contigs"])


################################################################################
## sORF prediction
################################################################################


rule filter_nt_contigs_to_short:
    input:
        all_contigs=ALL_CONTIGS,
        short_contigs=SHORT_CONTIGS,
    output:
        contigs300=temp(OUTPUT_DIR / "outputs/sORF/short_contigs/contigs300.fa"),
        all_short_contigs=OUTPUT_DIR / "sORF/short_contigs/short_contigs.fa",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit seq --max-len 300 -o {output.contigs300} {input.all_contigs}
        cat {input.short_contigs} {output.contigs300} > {output.all_short_contigs}
        """


# TER TODO: Add a rule for sORF prediction, either once smallesm is developed,
#           when there is an accurate sORF rnasamba model,
#           or using another tool from Singh & Roy.


rule filter_nt_contigs_to_long:
    input:
        all_contigs=ALL_CONTIGS,
    output:
        long_contigs=temp(OUTPUT_DIR / "sORF/long_contigs/contigs300.fa"),
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit seq --min-len 301 -o {output.long_contigs} {input.all_contigs}
        """


rule get_coding_contig_names:
    """
    Extract amino acid contig names and remove everything after the first period, 
    which are isoform labels.
    This file will be used to select all contigs that DO NOT encode ORFs, 
    according to transdecoder.
    """
    input:
        ORFS_AMINO_ACIDS,
    output:
        names=OUTPUT_DIR / "sORF/long_contigs/orfs_amino_acid_names.txt",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit seq -n {input} | sed 's/[.].*$//' > {output}
        """


rule filter_long_contigs_to_no_predicted_ORF:
    """
    Many of the contigs in the full transcriptome have predicted ORFs.
    The names of these contigs are recorded in the transdecoder input files (*pep & *cds, orfs_*).
    By definition, these contigs are not noncoding RNAs, 
    so they don't need to be considered for classification as long noncoding RNAs (lncRNA).
    This step removes the contigs that contain ORFs.
    """
    input:
        fa=rules.filter_nt_contigs_to_long.output.long_contigs,
        names=rules.get_coding_contig_names.output.names,
    output:
        fa=OUTPUT_DIR / "sORF/long_contigs/long_contigs_no_predicted_orf.fa",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit grep -v -f {input.names} {input.fa} -o {output.fa}
        """


rule download_rnasamba_model:
    """
    Place holder rule.
    For now, the workflow uses the model output by build_rnasamba_euk_model.snakefile, 
    which is available locally from running it.
    """
    output:
        model=OUTPUT_DIR / "models/rnasamba/build/3_model/eu_rnasamba.hdf5",
    shell:
        """
        curl -JLo {output.model} # TODO add URL for download
        """


rule pip_install_rnasamba_no_deps:
    """
    To take advantage of nvidia GPU on AWS instance 
    ("Deep Learning Base OSS Nvidia Driver GPU AMI (Ubuntu 20.04) 20240122" ami-07eb000b3340966b0),
    we need to install specific versions of tensorflow and other dependencies.
    This is accomplished in part in the envs/rnasamba.yml file, 
    however rnasamba itself is not installed there because we need to use the command:
    pip install --no-deps rnasamba
    and there is no way to specify the "--no-deps" flag in a yaml file.
    This rule installs rnasamba into the conda-generated environment.
    Note the output path is the same as that used by curate_datasets_and_build_models.snakefile, 
    as this conda env will already be present and configured if that snakefile has been executed.
    """
    output:
        pip="outputs/models/build/rnasamba/rnasamba_installed.txt",
    conda:
        "envs/rnasamba.yml"
    shell:
        """
        pip install 'nvidia-tensorflow~=1.15'
        pip install --no-deps rnasamba # used version 0.2.5
        touch {output}
        """


rule rnasamba:
    """
    The eukaryote_rnasamba.hdf5 model is only accurate on longer contigs.
    It assesses whether they are long noncoding RNAs.
    However, lncRNAs often have sORFs that encode peptides.
    This rule runs RNAsamba on longer contigs (>300nt) that were not predicted by transdecoder to 
    contain ORFs.
    """
    input:
        # TER TODO: update path when model is downloaded
        pip=rules.pip_install_rnasamba_no_deps.output.pip,
        model=rules.download_rnasamba_model.output.model,
        contigs=rules.filter_long_contigs_to_no_predicted_ORF.output.fa,
    output:
        tsv=OUTPUT_DIR / "sORF/long_contigs/rnasamba/classification.tsv",
        fa=OUTPUT_DIR / "sORF/long_contigs/rnasamba/predicted_proteins.fa",
    conda:
        "envs/rnasamba.yml"
    shell:
        """
        rnasamba classify -p {output.fa} {output.tsv} {input.contigs} {input.model}
        """


## TER TODO: predict sORFs from lncRNAs


################################################################################
## cleavage prediction
################################################################################


rule remove_stop_codon_asterisk_from_transdecoder_ORFs:
    input:
        ORFS_AMINO_ACIDS,
    output:
        faa=OUTPUT_DIR / "cleavage/preprocessing/noasterisk.faa",
    shell:
        """
        sed '/^[^>]/s/\*//g' {input} > {output}
        """


# Ribosomally synthesized and post-translationally modified peptide prediction


rule download_nlpprecursor_models:
    output:
        tar=INPUT_DIR / "models/nlpprecursor/nlpprecursor_models.tar.gz",
        model=INPUT_DIR / "models/nlpprecursor/models/annotation/model.p",
    params:
        outdir=INPUT_DIR / "models/nlpprecursor",
    shell:
        """
        curl -JLo {output.tar} https://github.com/magarveylab/NLPPrecursor/releases/download/1.0/nlpprecursor_models.tar.gz
        tar xf {output.tar} -C {params.outdir} 
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
        faa=rules.remove_stop_codon_asterisk_from_transdecoder_ORFs.output.faa,
        model=rules.download_nlpprecursor_models.output.model,
    output:
        tsv=OUTPUT_DIR / "cleavage/nlpprecursor/nlpprecursor_predictions.tsv",
        peptide=OUTPUT_DIR / "cleavage/nlpprecursor/nlpprecursor_peptides.fasta",
    params:
        modelsdir=INPUT_DIR / "models/nlpprecursor/models/",
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
        json=OUTPUT_DIR / "cleavage/deeppeptide/peptide_predictions.json",
    conda:
        "envs/deeppeptide.yml"
    params:
        outdir=OUTPUT_DIR / "cleavage/deeppeptide",
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
        propeptide=OUTPUT_DIR / "cleavage/deeppeptide/propeptides.faa",
        peptide=OUTPUT_DIR / "cleavage/deeppeptide/peptides.faa",
        tsv=OUTPUT_DIR / "cleavage/deeppeptide/predictions.tsv",
    conda:
        "envs/biopython.yml"
    shell:
        """
        python scripts/extract_deeppeptide_sequences.py {input.json} {input.faa} {output.propeptide} {output.peptide} {output.tsv}
        """


################################################################################
## Non-ribosomal peptide synthetase annotation
################################################################################


rule nrps_hmmsearch:
    """
    Uses hidden markov models to search for domains that are commonly annotated in NRPS genes.
    The input.hmm file is produced by download_nrps_hmm_profiles.snakefile and included in this repo. 
    """
    input:
        faa=rules.remove_stop_codon_asterisk_from_transdecoder_ORFs.output.faa,
        hmm=INPUT_DIR / "models/nrps/nrps.hmm",
    output:
        txt=OUTPUT_DIR / "nrps/hmmsearch/hmmsearch.txt",
        tbltsv=OUTPUT_DIR / "nrps/hmmsearch/hmmsearch.tbltsv",
        domtsv=OUTPUT_DIR / "nrps/hmmsearch/hmmsearch.domtsv",
    conda:
        "envs/hmmer.yml"
    threads: 4
    shell:
        """
        hmmsearch -o {output.txt} --tblout {output.tbltsv} --domtblout {output.domtsv} --cpu {threads} {input.hmm} {input.faa} 
        """


################################################################################
## Combine peptide predictions
################################################################################


# TER TODO: figure out if the NRPS results are significant and should be parsed and included
# TER TODO: figure out if nlpprecursor results need to be filtered
# TER TODO: add sORF predictions


rule combine_peptide_predictions:
    input:
        nlpprecursor=rules.nlpprecursor.output.peptide,
        deeppeptide=rules.extract_deeppeptide_sequences.output.peptide,
    output:
        peptide=OUTPUT_DIR / "annotation/combined_peptide_predictions/peptides.faa",
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
        db=INPUT_DIR / "databases/peptipedia.fasta.gz",
    shell:
        """
        curl -JLo {output} https://osf.io/dzycu/download 
        """


rule make_diamond_db_from_peptipedia_database:
    input:
        db=rules.download_peptipedia_database.output.db,
    output:
        db=OUTPUT_DIR / "annotation/peptipedia/0_diamond_db/peptipedia.dmnd",
    params:
        dbprefix=OUTPUT_DIR / "annotation/peptipedia/0_diamond_db/peptipedia",
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
        tsv=OUTPUT_DIR / "annotation/peptipedia/1_blastp/matches.tsv",
    params:
        dbprefix=OUTPUT_DIR / "annotation/peptipedia/0_diamond_db/peptipedia",
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
        tsv=OUTPUT_DIR / "annotation/deepsig/deepsig.tsv",
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
        tsv=OUTPUT_DIR / "annotation/characteristics/peptide_characteristics.tsv",
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
        # TER TODO: the authors of autopeptideml sent me these models.
        # They said they're working on uploading them.
        # Once they're available, I need to add a rule to download them and update the input here
        # to be the rules syntax
        model=INPUT_DIR
        / "models/autopeptideml/HPO_NegSearch_HP/{autopeptideml_model_name}_1/apml_config.json",
    output:
        tsv=OUTPUT_DIR / "annotation/autopeptideml/autopeptideml_{autopeptideml_model_name}.tsv",
    params:
        modelsdir=INPUT_DIR / "models/autopeptideml/HPO_NegSearch_HP/",
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
        autopeptideml=expand(
            rules.run_autopeptideml.output.tsv, autopeptideml_model_name=AUTOPEPTIDEML_MODEL_NAMES
        ),
        deepsig=rules.run_deepsig.output.tsv,
        peptipedia=rules.diamond_blastp_peptide_predictions_against_peptipedia_database.output.tsv,
        characteristics=rules.characterize_peptides.output.tsv,
    output:
        tsv=OUTPUT_DIR / "annotation/peptide_annotations.tsv",
    params:
        autopeptidemldir=OUTPUT_DIR / "annotation/autopeptideml/",
    conda:
        "envs/tidyverse.yml"
    shell:
        """
        Rscript scripts/combine_peptide_annotations.R \
            --nlpprecursor_path {input.nlpprecursor} \
            --deeppeptide_path {input.deeppeptide} \
            --autopeptideml_dir {params.autopeptidemldir} \
            --deepsig_path {input.deepsig} \
            --peptipedia_path {input.peptipedia} \
            --characteristics_path {input.characteristics} \
            --output_path {output.tsv}
        """


################################################################################
## Target rule all
################################################################################


rule all:
    default_target: True
    input:
        rules.rnasamba.output.tsv,
        rules.nlpprecursor.output.peptide,
        rules.extract_deeppeptide_sequences.output.peptide,
        rules.nrps_hmmsearch.output.tbltsv,


rule predict_sORF:
    """
    Defines a target rule for sORF prediction so a user can run only sORF prediction.
    snakemake predict_sORF --software-deployment-method conda -j 8 
    """
    input:
        rules.rnasamba.output.tsv,


rule predict_cleavage:
    """
    Defines a target rule for cleavage prediction so a user can run only cleavage prediction.
    snakemake predict_cleavage --software-deployment-method conda -j 8 
    """
    input:
        rules.combine_peptide_annotations.output.tsv,


rule predict_nrps:
    """
    Defines a target rule for nonribosomal peptide synthetase prediction so a user can run only NRPS prediction if they desire.
    snakemake predict_nrps --software-deployment-method conda -j 8 
    """
    input:
        rules.nrps_hmmsearch.output.tbltsv,
