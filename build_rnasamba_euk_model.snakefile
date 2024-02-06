import pandas as pd
import os
import re

metadata = pd.read_csv("inputs/models/rnasamba/build/train_data_links.tsv", sep = "\t")
GENOMES = metadata['genome'].unique().tolist()
RNA_TYPES = ['cdna', 'ncrna'] # inherits names from ensembl
VALIDATION_TYPES = ['mRNAs', 'ncRNAs'] # inherits names from https://github.com/cbl-nabi/RNAChallenge
SET_TYPES = ['coding', 'noncoding']
SET_NAMES = ['train', 'test', 'validation']

rule all:
    input: 
        expand("outputs/models/rnasamba/build/2_sequence_sets/{set_type}_{set_name}.fa", set_type = SET_TYPES, set_name = SET_NAMES)

rule download_ensembl_data:
    """
    Download ensembl cDNA and ncRNA files. 
    Ensembl annotates protein coding and non-coding RNA transcripts in their files.
    This information will be used to separate protein coding from non-coding RNAs to build an RNAsamba model.
    Note this download renames genome files from their names on ensembl to make them simpler to point to.
    See example transformations below:
    - cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -> cdna/Homo_sapiens.GRCh38.cdna.fa.gz (dropped "all.")
    - ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz   -> ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz (no change)
    """
    output: "inputs/ensembl/{rna_type}/{genome}.{rna_type}.fa.gz"
    run:
        genome_df = metadata.loc[(metadata['genome'] == wildcards.genome)]
        root_url = genome_df['root_url'].values[0]
        if wildcards.rna_type == "cdna":
            suffix = genome_df['cdna_suffix'].values[0]
        else:
            suffix = genome_df['ncrna_suffix'].values[0]
        
        url = root_url + suffix
        shell("curl -JLo {output} {url}")


rule extract_protein_coding_orfs_from_cdna:
    """
    Ensembl cDNA files consist of transcript sequences for actual and possible genes, including pseudogenes, NMD and the like.
    Transcripts in the cDNA files have headers like: >TRANSCRIPT_ID SEQTYPE LOCATION GENE_ID GENE_BIOTYPE TRANSCRIPT_BIOTYPE, where the gene_biotype and transcript_biotype both contain information about whether the gene is coding or not.
    """
    input: "inputs/ensembl/cdna/{genome}.cdna.fa.gz"
    output: "outputs/models/rnasamba/build/0_coding/{genome}.fa.gz"
    conda: "envs/seqkit.yml"
    shell:'''
    seqkit grep --use-regexp --by-name --pattern "transcript_biotype:protein_coding" -o {output} {input}
    '''

rule download_validation_data:
    output: "inputs/validation/rnachallenge/{validation_type}.fa.gz",
    shell:'''
    curl -JL https://raw.githubusercontent.com/cbl-nabi/RNAChallenge/main/RNAchallenge/{wildcards.validation_type}.fa | gzip > {output}
    '''

rule combine_sequences:
    input:
       coding = expand("outputs/models/rnasamba/build/0_coding/{genome}.fa.gz", genome = GENOMES),
       noncoding = expand("inputs/ensembl/ncrna/{genome}.ncrna.fa.gz", genome = GENOMES),
       validation = expand("inputs/validation/rnachallenge/{validation_type}.fa.gz", validation_type = VALIDATION_TYPES)
    output: "outputs/models/rnasamba/build/1_homology_reduction/all_sequences.fa.gz"
    shell:'''
    cat {input} > {output}
    '''
    
rule reduce_sequence_homology:
    """
    To reduce pollution between training and testing set, cluster sequences at 80% sequence identity.
    """
    input: "outputs/models/rnasamba/build/1_homology_reduction/all_sequences.fa.gz"
    output: 
        "outputs/models/rnasamba/build/1_homology_reduction/clustered_sequences_rep_seq.fasta",
        "outputs/models/rnasamba/build/1_homology_reduction/clustered_sequences_cluster.tsv"
    params: prefix = "outputs/models/rnasamba/build/1_homology_reduction/clustered_sequences"
    conda: "envs/mmseqs2.yml"
    shell:'''
    mmseqs easy-cluster {input} {params.prefix} tmp_mmseqs2 --min-seq-id 0.8 --cov-mode 1 --cluster-mode 2
    '''

rule grab_validation_set_names_and_lengths:
    input: "inputs/validation/rnachallenge/{validation_type}.fa.gz",
    output: 
        validation = "inputs/validation/rnachallenge/{validation_type}.fa",
        validation_fai = "inputs/validation/rnachallenge/{validation_type}.fa.fai",
    conda: "envs/seqkit.yml"
    shell:'''
    cat {input} | gunzip > {output.validation}
    seqkit faidx {output.validation}
    '''

# TER TODO: consider changing the ncrna input path to remove some duplication in the rules
rule grab_traintest_coding_names_and_lengths:
    input: expand("outputs/models/rnasamba/build/0_coding/{genome}.fa.gz", genome = GENOMES),
    output:
        coding = "outputs/models/rnasamba/build/2_sequence_sets/traintest/all_coding.fa",
        coding_fai = "outputs/models/rnasamba/build/2_sequence_sets/traintest/all_coding.fa.fai"
    conda: "envs/seqkit.yml"
    shell:'''
    cat {input} | gunzip > {output.coding}
    seqkit faidx {output.coding}
    '''

rule grab_traintest_noncoding_names_and_lengths:
    input: expand("inputs/ensembl/ncrna/{genome}.ncrna.fa.gz", genome = GENOMES),
    output:
        noncoding = "outputs/models/rnasamba/build/2_sequence_sets/traintest/all_noncoding.fa",
        noncoding_fai = "outputs/models/rnasamba/build/2_sequence_sets/traintest/all_noncoding.fa.fai"
    conda: "envs/seqkit.yml"
    shell:'''
    cat {input} | gunzip > {output.noncoding}
    seqkit faidx {output.noncoding}
    '''

rule process_sequences_into_nonoverlapping_sets:
    input: 
        traintest_fai = expand("outputs/models/rnasamba/build/2_sequence_sets/traintest/all_{set_type}.fa.fai", set_type = SET_TYPES), 
        validation_fai = expand("inputs/validation/rnachallenge/{validation_type}.fa.fai", validation_type = VALIDATION_TYPES),
        clusters = "outputs/models/rnasamba/build/1_homology_reduction/clustered_sequences_cluster.tsv"
    output:
        # TER TODO: the set_name/set_types could be born here, but I need to check the order in which they would be built to make sure files aren't named the wrong thing. actually i could do this by number of lines in the final file. Also need to update the R script to point to the right thing
        coding_train = "outputs/models/rnasamba/build/2_sequence_sets/coding_train.txt",
        coding_test = "outputs/models/rnasamba/build/2_sequence_sets/coding_test.txt",
        noncoding_train = "outputs/models/rnasamba/build/2_sequence_sets/noncoding_train.txt",
        noncoding_test = "outputs/models/rnasamba/build/2_sequence_sets/noncoding_test.txt",
        coding_validation = "outputs/models/rnasamba/build/2_sequence_sets/coding_validation.txt",
        noncoding_validation = "outputs/models/rnasamba/build/2_sequence_sets/noncoding_validation.txt"
    conda: "envs/tidyverse.yml"
    script: "scripts/process_sequences_into_nonoverlapping_sets.R"

rule filter_sequence_sets:
    input:
        fa = "outputs/models/rnasamba/build/1_homology_reduction/clustered_sequences_rep_seq.fasta",
        names = "outputs/models/rnasamba/build/2_sequence_sets/{set_type}_{set_name}.txt"
    output: "outputs/models/rnasamba/build/2_sequence_sets/{set_type}_{set_name}.fa"
    conda: "envs/seqtk.yml"
    shell:'''
    seqtk subseq {input.fa} {input.names} > {output}
    '''

##################################################################
## Get sequence statistics
##################################################################

rule get_sequence_statistics:
    input: "outputs/models/rnasamba/build/2_sequence_sets/{set_type}_{set_name}.fa"
    output: "outputs/models/rnasamba/build/stats/full/{set_type}_{set_name}.tsv"
    conda: "envs/seqkit.yml"
    shell:'''
    seqkit stats --all -o {output} -T {input}
    '''

rule get_sequence_statistics_less_than_300_nt:
    input: "outputs/models/rnasamba/build/2_sequence_sets/{set_type}_{set_name}.fa"
    output: "outputs/models/rnasamba/build/stats/300nt_or_less/{set_type}_{set_name}.tsv"
    conda: "envs/seqkit.yml"
    shell:'''
    seqkit seq --max-len 300 {input} | seqkit stats --all -T -o {output}
    '''

def read_tsv_files(file_paths):
    all_data = [] 
    
    for file_path in file_paths:
        file_name = re.sub("outputs/models/rnasamba/build/stats/", "", file_path)
        file_type = re.sub("/.*", "", file_name)
        if "cdna" in file_name:
            rna_type = "cnda"
        elif "ncrna" in file_name:
            rna_type = "ncrna"
        else:
            rna_type = "protein_coding"
        rm_from_genome_str1 = "." + rna_type + ".tsv" 
        genome = re.sub(rm_from_genome_str1, "", file_name)
        rm_from_genome_str2 = file_type + "/"
        genome = re.sub(rm_from_genome_str2, "", genome)
        df = pd.read_csv(file_path, sep='\t', header=0)
        df['file_type'] = file_type
        df['rna_type'] = rna_type
        df['genome'] = genome 
        all_data.append(df)
    
    # Concatenate all dataframes
    combined_df = pd.concat(all_data, ignore_index=True)
    return combined_df


rule combine_stats_files:
    input:
        expand("outputs/models/rnasamba/build/stats/full/{genome}.{rna_type}.tsv", rna_type = RNA_TYPES, genome = GENOMES),
        expand("outputs/models/rnasamba/build/stats/300nt_or_less/{genome}.ncrna.tsv", genome = GENOMES),
        expand("outputs/models/rnasamba/build/stats/300nt_or_less/{genome}.protein_coding.tsv", genome = GENOMES),
        expand("outputs/models/rnasamba/build/stats/protein_coding/{genome}.protein_coding.tsv", genome = GENOMES)
    output: "outputs/models/rnasamba/build/stats.tsv"
    run:
        combined_df = read_tsv_files(input)
        combined_df.to_csv(str(output), sep = "\t")
