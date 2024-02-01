################################################################################
## Input file descriptions
################################################################################

# All input files are produced by reads2transcriptome.
# While this is not a strict requirement, using these output files gives us the ability to look at very short contiguous sequences.
# TER TODO: update names of output files based on Nextflow of reads2transcriptome.
# - small_contigs.fa: contigs that are shorter than X nucleotides, which do not progress through the assembly pipeline. These may contain sORFs and so are included.
# - orthofuser_final_clean.fa.transdecoder.pep: predicted ORFs translated into amino acids. Output by dammit with ORF prediction by transdecoder. Used for cleavage peptide prediction and annotation of nonribosomal peptide synthetases.
# - orthofuser_final_clean.fa.transdecoder.cds: predicted ORFs as nucleotide sequences. Output by dammit with ORF prediction by transdecoder. Used to compare peptide nucleotide sequences (clustering, dn/ds estimation, etc.).
# - orthofuser_final_clean.fa.dammit.fasta: all contigs as nucleotide sequences. Used to identify contigs shorter than X nucleotides (300) to scan for sORFs and to predict lncRNAs, which may have sORFs embedded in them.
 
short_contigs = config.get("short_contigs", "inputs/reads2transcriptome_outputs/small_contigs.fa")
orfs_amino_acids = config.get("orfs_amino_acids", "inputs/reads2transcriptome_outputs/orthofuser_final_clean.fa.transdecoder.pep")
orfs_nucleotides = config.get("orfs_nucleotides", "inputs/reads2transcriptome_outputs/orthofuser_final_clean.fa.transdecoder.cds")
all_contigs = config.get("all_contigs", "inputs/reads2transcriptome_outputs/orthofuser_final_clean.fa.dammit.fasta")

rule all:
    input:

################################################################################
## sORF prediction
################################################################################

LENGTH = ['short', 'long']

rule sORF:
    """
    Defines a target rule for sORF prediction so a user can run only sORF prediction if they desire.
    snakemake sORF --software-deployment-method conda -j 8 
    """
    input:

rule filter_nt_contigs_to_short:
    input:
        all_contigs = all_contigs,
        short_contigs = short_contigs
    output:
        contigs300 = temp("outputs/sORF/short_contigs/contigs300.fa"),
        all_short_contigs = "outputs/sORF/short_contigs/short_contigs.fa"
    conda: "envs/seqkit.yml"
    shell:'''
    seqkit seq --max-len 300 -o {output.contigs300} {input.all_contigs}
    cat {input.short_contigs} {output.contigs300} > {output.all_short_contigs}
    '''

rule filter_nt_contigs_to_long:
    input: all_contigs = all_contigs
    # TER todo make temp()
    output: "outputs/sORF/long_contigs/contigs300.fa"
    conda: "envs/seqkit.yml"
    shell:'''
    seqkit seq --min-len 301 -o {output} {input.all_contigs}
    '''

rule get_coding_contig_names:
    """
    Extract amino acid contig names and remove everything after the first period, which are isoform labels.
    """
    input: orfs_amino_acids,
    output: "outputs/sORF/long_contigs/orfs_amino_acid_names.txt"
    conda: "envs/seqkit.yml"
    shell:'''
    seqkit seq -n {input} | sed 's/[.].*$//' > {output}
    '''

rule filter_long_contigs_to_noncoding:
    '''
    Many of the contigs in the full transcriptome have predicted ORFs.
    The names of these contigs are recorded in the transdecoder input files (*pep and *cds, orfs_*).
    By definition, these contigs are not noncoding RNAs, so they don't need to be considered for classification as long noncoding RNAs (lncRNA).
    This step removes the contigs that contain ORFs.
    '''
    input: 
        fa = "outputs/sORF/long_contigs/contigs300.fa",
        names = "outputs/sORF/long_contigs/orfs_amino_acid_names.txt"
    output: "outputs/sORF/long_contigs/long_contigs.fa"
    conda: "envs/seqkit.yml"
    shell:'''
    seqkit grep -v -f {input.names} {input.fa} -o {output}
    '''

rule download_rnasamba_model:
    output: "inputs/models/rnasamba/full_length_weights.hdf"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://raw.githubusercontent.com/apcamargo/RNAsamba/master/data/full_length_weights.hdf
    '''

rule rnasamba:
    """
    RNAsamba will run on both the short and long contigs.
    For short contigs, it assess their coding potential.
    Short contigs that have coding potential may encode peptides.
    For long contigs, it assesses whether they are long noncoding RNAs.
    LncRNAs often have sORFs that encode peptides.
    """
    input:
        model = "inputs/models/rnasamba/full_length_weights.hdf",
        contigs = "outputs/sORF/{length}_contigs/{length}_contigs.fa"
    output:
        tsv = "outputs/sORF/{length}_contigs/rnasamba/classification.tsv",
        fa  = "outputs/sORF/{length}_contigs/rnasamba/predicted_proteins.fa"
    conda: "envs/rnasamba.yml"
    shell:'''
    rnasamba classify -p {output.fa} {output.tsv} {input.all_contigs} {input.model}
    '''

