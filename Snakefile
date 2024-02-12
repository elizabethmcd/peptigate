import os
from pathlib import Path

################################################################################
## Configuration
################################################################################

# Default pipeline configuration parameters are in the config file.
# If you create a new yml file and use the --configfile flag, options in that new file overwrite the defaults.
configfile: "./config.yml"

input_dir = Path(config["input_dir"])
output_dir = Path(config["output_dir"])

short_contigs = Path(config["short_contigs"])
orfs_amino_acids = Path(config["orfs_amino_acids"])
orfs_nucleotides = Path(config["orfs_nucleotides"])
all_contigs = Path(config["all_contigs"])



################################################################################
## sORF prediction
################################################################################


rule sORF:
    """
    Defines a target rule for sORF prediction so a user can run only sORF prediction if they desire.
    snakemake sORF --software-deployment-method conda -j 8 
    """
    input:
        rnasamba="outputs/sORF/long_contigs/rnasamba/classification.tsv",


rule filter_nt_contigs_to_short:
    input:
        all_contigs=all_contigs,
        short_contigs=short_contigs,
    output:
        contigs300=temp("outputs/sORF/short_contigs/contigs300.fa"),
        all_short_contigs="outputs/sORF/short_contigs/short_contigs.fa",
    conda:
        "envs/seqkit.yml"
    shell:
        """
    seqkit seq --max-len 300 -o {output.contigs300} {input.all_contigs}
    cat {input.short_contigs} {output.contigs300} > {output.all_short_contigs}
    """


# TER TODO: Add a rule for sORF prediction, either once smallesm is developed, when there is an accurate sORF rnasamba model, or using another tool from Singh & Roy.


rule filter_nt_contigs_to_long:
    input:
        all_contigs=all_contigs,
    output:
        long_contigs=temp("outputs/sORF/long_contigs/contigs300.fa"),
    conda:
        "envs/seqkit.yml"
    shell:
        """
    seqkit seq --min-len 301 -o {output} {input.all_contigs}
    """


rule get_coding_contig_names:
    """
    Extract amino acid contig names and remove everything after the first period, which are isoform labels.
    This file will be used to select all contigs that DO NOT encode ORFs, according to transdecoder.
    """
    input:
        orfs_amino_acids,
    output:
        names = "outputs/sORF/long_contigs/orfs_amino_acid_names.txt",
    conda:
        "envs/seqkit.yml"
    shell:
        """
    seqkit seq -n {input} | sed 's/[.].*$//' > {output}
    """


rule filter_long_contigs_to_no_predicted_ORF:
    """
    Many of the contigs in the full transcriptome have predicted ORFs.
    The names of these contigs are recorded in the transdecoder input files (*pep and *cds, orfs_*).
    By definition, these contigs are not noncoding RNAs, so they don't need to be considered for classification as long noncoding RNAs (lncRNA).
    This step removes the contigs that contain ORFs.
    """
    input:
        fa = rules.filter_nt_contigs_to_long.output.long_contigs, 
        names = rules.get_coding_contig_names.output.names
    output:
        fa = "outputs/sORF/long_contigs/long_contigs_no_predicted_orf.fa",
    conda:
        "envs/seqkit.yml"
    shell:
        """
    seqkit grep -v -f {input.names} {input.fa} -o {output.fa}
    """


rule download_rnasamba_model:
    """
    Place holder rule.
    For now, the workflow uses the model output by build_rnasamba_euk_model.snakefile, which is available locally from running it.
    """
    output:
        model = "outputs/models/rnasamba/build/3_model/eu_rnasamba.hdf5",
    shell:
        """
    curl -JLo {output.model} # TODO add URL for download
    """


rule rnasamba:
    """
    The eu_rnasamba.hdf5 model is only accurate on longer contigs.
    It assesses whether they are long noncoding RNAs.
    However, lncRNAs often have sORFs that encode peptides.
    This rule runs RNAsamba on longer contigs (>300nt) that were not predicted by transdecoder to contain ORFs.
    """
    input:
        # TER TODO: update path when model is downloaded
        model=rules.download_rnasamba_model.output.model,
        contigs=rules.filter_long_contigs_to_no_predicted_ORF.output.fa
    output:
        tsv="outputs/sORF/long_contigs/rnasamba/classification.tsv",
        fa="outputs/sORF/long_contigs/rnasamba/predicted_proteins.fa",
    conda:
        "envs/rnasamba.yml"
    shell:
        """
    rnasamba classify -p {output.fa} {output.tsv} {input.contigs} {input.model}
    """


## TER TODO: predict sORFs from lncRNAs


################################################################################
## Target rule all 
################################################################################

rule all:
    default_target: True
    input:
        rules.rnasamba.output.tsv
