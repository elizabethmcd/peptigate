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
        "output.txt",

################################################################################
## sORF prediction
################################################################################

rule sORF:
    """
    Defines a target rule for sORF prediction so a user can run only sORF prediction if they desire.
    snakemake sORF --software-deployment-method conda -j 8 
    """
    input:

rule filter_nt_contigs_by_length:
    output:
        "output.txt",
    shell:'''
    
    '''
