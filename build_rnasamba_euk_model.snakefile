import pandas as pd
import os
import re

metadata = pd.read_csv("inputs/models/rnasamba/build/train_data_links.tsv", sep = "\t")
GENOMES = metadata['genome'].unique().tolist()
RNA_TYPES = ['cdna', 'ncrna'] # inherits names from ensembl
VALIDATION_TYPES = ['mRNAs', 'ncRNAs'] # inherits names from https://github.com/cbl-nabi/RNAChallenge

rule all:
    input: 
        "outputs/models/rnasamba/build/stats.tsv",
        "outputs/models/rnasamba/build/2_sequence_sets/protein_coding_traintest.fa",
        "outputs/models/rnasamba/build/2_sequence_sets/non_coding_traintest.fa",
        expand("outputs/models/rnasamba/build/2_sequence_sets/validation/{validation_type}.fa", validation_type = VALIDATION_TYPES)

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
    output: "outputs/models/rnasamba/build/0_protein_coding/{genome}.fa.gz"
    conda: "envs/seqkit.yml"
    shell:'''
    seqkit grep --use-regexp --by-name --pattern "transcript_biotype:protein_coding" -o {output} {input}
    '''

rule download_validation_data:
    output: "inputs/validation/rnachallenge/{validation_type}.fa",
    shell:'''
    curl -JLo {output} https://raw.githubusercontent.com/cbl-nabi/RNAChallenge/main/RNAchallenge/{wildcards.validation_type}.fa
    '''

rule combine_sequences:
    input:
       coding = expand("outputs/models/rnasamba/build/0_protein_coding/{genome}.fa.gz", genome = GENOMES),
       noncoding = expand("inputs/ensembl/ncrna/{genome}.ncrna.fa.gz", genome = GENOMES),
       validation = expand("inputs/validation/rnachallenge/{validation_type}.fa", validation_type = VALIDATION_TYPES)
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

rule grab_representative_sequence_names:
    input: "outputs/models/rnasamba/build/1_homology_reduction/clustered_sequences_rep_seq.fasta",
    output: "outputs/models/rnasamba/build/1_homology_reduction/clustered_sequences_rep_seq_names.txt"
    shell:'''
    awk 'sub(/^>/, "")' {input} > {output}
    '''

rule filter_validation:
    input:
        fa = "inputs/validation/rnachallenge/{validation_type}.fa",
        names = "outputs/models/rnasamba/build/1_homology_reduction/clustered_sequences_rep_seq_names.txt"
    output: "outputs/models/rnasamba/build/2_sequence_sets/validation/{validation_type}.fa"
    conda: "envs/seqtk.yml"
    shell:'''
    seqtk subseq {input.fa} {input.names} > {output}
    '''

rule filter_traintest_coding:
    input:
        protein_coding = expand("outputs/models/rnasamba/build/0_protein_coding/{genome}.fa.gz", genome = GENOMES),
        names = "outputs/models/rnasamba/build/1_homology_reduction/clustered_sequences_rep_seq_names.txt"
    output:
        protein_coding = temp("outputs/models/rnasamba/build/2_sequence_sets/all_protein_coding.fa.gz"),
        protein_coding_filt = "outputs/models/rnasamba/build/2_sequence_sets/protein_coding_traintest.fa"
    conda: "envs/seqtk.yml"
    shell:'''
    cat {input.protein_coding} > {output.protein_coding}
    seqtk subseq {output.protein_coding} {input.names} > {output.protein_coding_filt}
    '''


rule filter_traintest_noncoding:
    input:
        non_coding = expand("inputs/ensembl/ncrna/{genome}.ncrna.fa.gz", genome = GENOMES),
        names = "outputs/models/rnasamba/build/1_homology_reduction/clustered_sequences_rep_seq_names.txt"
    output:
        non_coding = temp("outputs/models/rnasamba/build/2_sequence_sets/all_non_coding.fa.gz"),
        non_coding_filt = "outputs/models/rnasamba/build/2_sequence_sets/non_coding_traintest.fa"
    conda: "envs/seqtk.yml"
    shell:'''
    cat {input.non_coding} > {output.non_coding}
    seqtk subseq {output.non_coding} {input.names} > {output.non_coding_filt}
    '''

##################################################################
## Get sequence statistics
##################################################################

rule get_sequence_statistics:
    input: "inputs/ensembl/{rna_type}/{genome}.{rna_type}.fa.gz"
    output: "outputs/models/rnasamba/build/stats/full/{genome}.{rna_type}.tsv"
    conda: "envs/seqkit.yml"
    shell:'''
    seqkit stats --all -o {output} -T {input}
    '''

rule get_sequence_statistics_less_than_300_nt_ncrna:
    input: "inputs/ensembl/ncrna/{genome}.ncrna.fa.gz"
    output: "outputs/models/rnasamba/build/stats/300nt_or_less/{genome}.ncrna.tsv"
    conda: "envs/seqkit.yml"
    shell:'''
    seqkit seq --max-len 300 {input} | seqkit stats --all -T -o {output}
    '''
rule get_sequence_statistics_less_than_300nt_protein_coding:
    input: "outputs/models/rnasamba/build/0_protein_coding/{genome}.fa.gz"
    output: "outputs/models/rnasamba/build/stats/300nt_or_less/{genome}.protein_coding.tsv"
    conda: "envs/seqkit.yml"
    shell:'''
    seqkit seq --max-len 300 {input} | seqkit stats --all -T -o {output}
    '''

rule get_sequence_statistics_protein_coding:
    input: "outputs/models/rnasamba/build/0_protein_coding/{genome}.fa.gz" 
    output: "outputs/models/rnasamba/build/stats/protein_coding/{genome}.protein_coding.tsv"
    conda: "envs/seqkit.yml"
    shell:'''
    seqkit stats --all -T -o {output} {input}
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
