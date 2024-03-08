import pandas as pd
import os
import re

metadata = pd.read_csv("inputs/models/datasets/train_data_links.tsv", sep="\t")
GENOMES = metadata["genome"].unique().tolist()
RNA_TYPES = ["cdna", "ncrna"]  # inherits names from ensembl
VALIDATION_TYPES = [
    "mRNAs",
    "ncRNAs",
]  # inherits names from https://github.com/cbl-nabi/RNAChallenge
# Note that for variables CODING_TYPES and DATASET_TYPES, the order matters;
# The snakemake expand() function maintains the order of these arguments.
# This behavior is used later in the snakefile to specify inputs and outputs of different scripts/shell commands.
# The order of CODING_TYPES must correspond to the order of the positional args that are later passed to rnasamba.
# In addition, the order of CODING_TYPES and DATASET_TYPES is used to declare outputs in the rule/R script process_sequences_into_nonoverlapping_sets.
CODING_TYPES = ["coding", "noncoding"]
DATASET_TYPES = ["train", "validation"]


rule all:
    input:
        "outputs/models/datasets/3_stats/set_summary.tsv",
        expand(
            "outputs/models/build/plmutils/3_predict/{coding_type}_validation_predictions.csv",
            coding_type=CODING_TYPES,
        ),


rule download_ensembl_data:
    """
    Download ensembl cDNA and ncRNA files.
    Ensembl annotates protein coding and non-coding RNA transcripts in their files.
    This information will be used to separate protein coding from non-coding RNAs.
    Note this download script changes the output file names to make them simpler to point to.
    See example transformations below:
    - cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -> cdna/Homo_sapiens.GRCh38.cdna.fa.gz (dropped "all.")
    - ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz   -> ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz (no change)
    """
    output:
        "inputs/models/datasets/ensembl/{rna_type}/{genome}.{rna_type}.fa.gz",
    run:
        genome_df = metadata.loc[(metadata["genome"] == wildcards.genome)]
        root_url = genome_df["root_url"].values[0]
        if wildcards.rna_type == "cdna":
            suffix = genome_df["cdna_suffix"].values[0]
        else:
            suffix = genome_df["ncrna_suffix"].values[0]

        url = root_url + suffix
        shell("curl -JLo {output} {url}")


rule extract_protein_coding_orfs_from_cdna:
    """
    Ensembl cDNA files consist of transcript sequences for actual and possible genes,
    including pseudogenes, NMD and the like.
    Transcripts in the cDNA files have headers like: 
    >TRANSCRIPT_ID SEQTYPE LOCATION GENE_ID GENE_BIOTYPE TRANSCRIPT_BIOTYPE,
    where gene_biotype and transcript_biotype contain information about whether the gene is coding.
    """
    input:
        "inputs/models/datasets/ensembl/cdna/{genome}.cdna.fa.gz",
    output:
        "outputs/models/datasets/0_coding/{genome}.fa.gz",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit grep \
            --use-regexp \
            --by-name \
            --pattern "transcript_biotype:protein_coding" \
            -o {output} \
            {input}
        """


rule download_validation_data:
    output:
        "inputs/models/datasets/validation/rnachallenge/{validation_type}.fa.gz",
    shell:
        """
        curl -JL https://raw.githubusercontent.com/cbl-nabi/RNAChallenge/main/RNAchallenge/{wildcards.validation_type}.fa | gzip > {output}
        """


rule combine_sequences:
    input:
        coding=expand("outputs/models/datasets/0_coding/{genome}.fa.gz", genome=GENOMES),
        noncoding=expand(
            "inputs/models/datasets/ensembl/ncrna/{genome}.ncrna.fa.gz", genome=GENOMES
        ),
        validation=expand(
            "inputs/models/datasets/validation/rnachallenge/{validation_type}.fa.gz",
            validation_type=VALIDATION_TYPES,
        ),
    output:
        "outputs/models/datasets/1_homology_reduction/all_sequences.fa",
    shell:
        """
        cat {input} | gunzip > {output}
        """


rule grab_all_sequence_names_and_lengths:
    input:
        "outputs/models/datasets/1_homology_reduction/all_sequences.fa",
    output:
        "outputs/models/datasets/1_homology_reduction/all_sequences.fa.seqkit.fai",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit faidx -f {input}
        """


rule reduce_sequence_homology:
    """
    Cluster sequences at 80% sequence identity to reduce pollution between train and test data sets.
    """
    input:
        "outputs/models/datasets/1_homology_reduction/all_sequences.fa",
    output:
        "outputs/models/datasets/1_homology_reduction/clustered_sequences_rep_seq.fasta",
        "outputs/models/datasets/1_homology_reduction/clustered_sequences_cluster.tsv",
    params:
        prefix="outputs/models/datasets/1_homology_reduction/clustered_sequences",
    conda:
        "envs/mmseqs2.yml"
    shell:
        """
        mmseqs easy-cluster {input} {params.prefix} tmp_mmseqs2 --min-seq-id 0.8 --cov-mode 1 --cluster-mode 2
        """


rule grab_validation_set_names_and_lengths:
    """
    The train/test data set sequences are identifiable by the genome information in the header,
    which is consistently formatted by Ensembl.
    The same is not true for the validation data.
    This rule grabs the validation sequence header names so they can be separated from the train/test sets.
    """
    input:
        "inputs/models/datasets/validation/rnachallenge/{validation_type}.fa.gz",
    output:
        validation="inputs/models/datasets/validation/rnachallenge/{validation_type}.fa",
        validation_fai="inputs/models/datasets/validation/rnachallenge/{validation_type}.fa.seqkit.fai",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        cat {input} | gunzip > {output.validation}
        seqkit faidx -f {output.validation}
        """


rule process_sequences_into_nonoverlapping_sets:
    input:
        metadata="inputs/models/datasets/train_data_links.tsv",
        all_fai="outputs/models/datasets/1_homology_reduction/all_sequences.fa.seqkit.fai",
        validation_fai=expand(
            "inputs/models/datasets/validation/rnachallenge/{validation_type}.fa.seqkit.fai",
            validation_type=VALIDATION_TYPES,
        ),
        clusters="outputs/models/datasets/1_homology_reduction/clustered_sequences_cluster.tsv",
    output:
        expand(
            "outputs/models/datasets/2_sequence_sets/{coding_type}_{dataset_type}.txt",
            coding_type=CODING_TYPES,
            dataset_type=DATASET_TYPES,
        ),
    conda:
        "envs/tidyverse.yml"
    script:
        "scripts/process_sequences_into_nonoverlapping_sets.R"


rule filter_sequence_sets:
    input:
        fa="outputs/models/datasets/1_homology_reduction/clustered_sequences_rep_seq.fasta",
        names="outputs/models/datasets/2_sequence_sets/{coding_type}_{dataset_type}.txt",
    output:
        "outputs/models/datasets/2_sequence_sets/{coding_type}_{dataset_type}.fa",
    conda:
        "envs/seqtk.yml"
    shell:
        """
        seqtk subseq {input.fa} {input.names} > {output}
        """


##################################################################
## Build plm-utils model
##################################################################


rule plmutils_translate:
    input:
        "outputs/models/datasets/2_sequence_sets/{coding_type}_{dataset_type}.fa",
    output:
        "outputs/models/build/plmutils/0_translate/{coding_type}_{dataset_type}.fa",
    conda:
        "envs/plmutils.yml"
    shell:
        """
        plmutils translate --longest-only --output-filepath {output}_tmp {input}
        # If this is the training data set, filter to sequences that are less than 100 amino acids.
        # Otherwise, keep all sequences.
        # We only want to train with short sequences since this our use case (detecting sORFs).
        # This is combine with this rule so that the output file naming scheme is consistent,
        # which allows us not to duplicate snakemake rules for downstream plmutils commands.
        if [ {wildcards.dataset_type} == "train" ]; then
            seqkit seq --max-len 100 -o {output} {output}_tmp
        else
            mv {output}_tmp {output}
        fi
        """


rule plmutils_embed:
    input:
        "outputs/models/build/plmutils/0_translate/{coding_type}_{dataset_type}.fa",
    output:
        "outputs/models/build/plmutils/1_embeddings/{coding_type}_{dataset_type}.npy",
    conda:
        "envs/plmutils.yml"
    shell:
        """
        plmutils embed --model-name esm2_t6_8M_UR50D \
            --layer-ind -1 \
            --output-filepath {output} \
            {input}
        """


rule plmutils_train:
    input:
        expand(
            "outputs/models/build/plmutils/1_embeddings/{coding_type}_train.npy",
            coding_type=CODING_TYPES,
        ),
    output:
        "outputs/models/build/plmutils/2_model/classifier.joblib",
    params:
        modeldir="outputs/models/build/plmutils/2_model/",
    conda:
        "envs/plmutils.yml"
    shell:
        """
        plmutils train --positive-class-filepath {input[0]} \
            --negative-class-filepath {input[1]} \
            --model-dirpath {params.modeldir}
        """


rule plmutils_predict_on_validation:
    input:
        embeddings="outputs/models/build/plmutils/1_embeddings/{coding_type}_validation.npy",
        fasta="outputs/models/build/plmutils/0_translate/{coding_type}_validation.fa",
        model="outputs/models/build/plmutils/2_model/classifier.joblib",
    output:
        "outputs/models/build/plmutils/3_predict/{coding_type}_validation_predictions.csv",
    params:
        modeldir="outputs/models/build/plmutils/2_model/",
    conda:
        "envs/plmutils.yml"
    shell:
        """
        plmutils predict --model-dirpath {params.modeldir} \
            --embeddings-filepath {input.embeddings} \
            --fasta-filepath {input.fasta} \
            --output-filepath {output}
        """


##################################################################
## Get sequence statistics
##################################################################


rule get_sequence_descriptors:
    input:
        "outputs/models/datasets/2_sequence_sets/{coding_type}_{dataset_type}.fa",
    output:
        "outputs/models/datasets/2_sequence_sets/{coding_type}_{dataset_type}.fa.seqkit.fai",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit faidx -f {input}
        """


rule calculate_sequence_statistics:
    input:
        expand(
            "outputs/models/datasets/2_sequence_sets/{coding_type}_{dataset_type}.fa.seqkit.fai",
            coding_type=CODING_TYPES,
            dataset_type=DATASET_TYPES,
        ),
    output:
        set_summary="outputs/models/datasets/3_stats/set_summary.tsv",
        set_length_summary="outputs/models/datasets/3_stats/set_length_summary.tsv",
        set_length_genome_summary="outputs/models/datasets/3_stats/set_length_genome_summary.tsv",
    conda:
        "envs/tidyverse.yml"
    script:
        "scripts/calculate_sequence_statistics.R"
