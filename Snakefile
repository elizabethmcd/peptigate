import os
from pathlib import Path

################################################################################
## Configuration
################################################################################


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
        tsv=OUTPUT_DIR / "cleavage/nlpprecursor/nlpprecursor_ripp_predictions.tsv",
    params:
        modelsdir=INPUT_DIR / "models/nlpprecursor/models/",
    conda:
        "envs/nlpprecursor.yml"
    shell:
        """
        python scripts/run_nlpprecursor.py {params.modelsdir} {input.faa} {output}
        """


# General Cleavage peptide prediction


rule clone_deeppeptide:
    output:
        src="cloned_repositories/DeepPeptide/LICENSE",
    shell:
        """
        cd cloned_repositories
        git clone https://github.com/fteufel/DeepPeptide.git
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
        outdir1=OUTPUT_DIR / "cleavage/deeppeptide/",
        outdir2=OUTPUT_DIR / "cleavage/",
    shell:
        """
        cd cloned_repositories/DeepPeptide/predictor && python3 predict.py --fastafile ../../../{input.faa} --output_dir {params.outdir1} --output_fmt json
        mv {params.outdir1} ../../../{params.outdir2}
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
    conda:
        "envs/deeppeptide.yml"
    shell:
        """
        python scripts/extract_deeppeptide_sequences.py {input.json} {input.faa} {output.propeptide} {output.peptide}
        """


################################################################################
## Target rule all
################################################################################


rule all:
    default_target: True
    input:
        rules.rnasamba.output.tsv,
        rules.nlpprecursor.output.tsv,
        rules.extract_deeppeptide_sequences.output.peptide,


rule sORF:
    """
    Defines a target rule for sORF prediction so a user can run only sORF prediction.
    snakemake sORF --software-deployment-method conda -j 8 
    """
    input:
        rules.rnasamba.output.tsv,


rule cleavage:
    """
    Defines a target rule for cleavage prediction so a user can run only cleavage prediction.
    snakemake cleavage --software-deployment-method conda -j 8 
    """
    input:
        rules.nlpprecursor.output.tsv,
        rules.extract_deeppeptide_sequences.output.peptide,
