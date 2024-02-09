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
    input: "outputs/sORF/long_contigs/rnasamba/classification.tsv",

################################################################################
## sORF prediction
################################################################################

rule sORF:
    """
    Defines a target rule for sORF prediction so a user can run only sORF prediction if they desire.
    snakemake sORF --software-deployment-method conda -j 8 
    """
    input: "outputs/sORF/long_contigs/rnasamba/classification.tsv",

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

# TER TODO: Add a rule for sORF prediction, either once smallesm is developed, when there is an accurate sORF rnasamba model, or using another tool from Singh & Roy.

rule filter_nt_contigs_to_long:
    input: all_contigs = all_contigs
    output: temp("outputs/sORF/long_contigs/contigs300.fa")
    conda: "envs/seqkit.yml"
    shell:'''
    seqkit seq --min-len 301 -o {output} {input.all_contigs}
    '''

rule get_coding_contig_names:
    """
    Extract amino acid contig names and remove everything after the first period, which are isoform labels.
    This file will be used to select all contigs that DO NOT encode ORFs, according to transdecoder.
    """
    input: orfs_amino_acids,
    output: "outputs/sORF/long_contigs/orfs_amino_acid_names.txt"
    conda: "envs/seqkit.yml"
    shell:'''
    seqkit seq -n {input} | sed 's/[.].*$//' > {output}
    '''

rule filter_long_contigs_to_no_predicted_ORF:
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

#rule download_rnasamba_model:
#    """
#    Place holder rule.
#    For now, the workflow uses the model output by build_rnasamba_euk_model.snakefile, which is available locally from running it. 
#    """
#    output: "inputs/models/rnasamba/eu_rnasamba.hdf5"
#    conda: "envs/wget.yml"
#    shell:'''
#    wget -O {output} 
#    '''

rule rnasamba:
    """
    The eu_rnasamba.hdf5 model is only accurate on longer contigs.
    It assesses whether they are long noncoding RNAs.
    However, lncRNAs often have sORFs that encode peptides.
    This rule runs RNAsamba on longer contigs (>300nt) that were not predicted by transdecoder to contain ORFs.
    """
    input:
        # TER TODO: update path when model is downloaded
        model = "outputs/models/rnasamba/build/3_model/eu_rnasamba.hdf5" 
        contigs = "outputs/sORF/long_contigs/long_contigs.fa"
    output:
        tsv = "outputs/sORF/long_contigs/rnasamba/classification.tsv",
        fa  = "outputs/sORF/long_contigs/rnasamba/predicted_proteins.fa"
    conda: "envs/rnasamba.yml"
    shell:'''
    rnasamba classify -p {output.fa} {output.tsv} {input.all_contigs} {input.model}
    '''

## TER TODO: predict sORFs from lncRNAs

