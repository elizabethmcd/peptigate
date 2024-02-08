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
orfs_amino_acids = config.get("orfs_amino_acids", "inputs/reads2transcriptome_outputs/orthofuser_final_clean.fa.transdecoder.pep_head")
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
    """
    The RNAsamba model was trained on human but had good predictive accuracy for distantly related model organisms.
    Accuracy did drop for drosophila and zebrafish, however the authors suggest that this might be due to inaccuracies in annotations in these two species.
    While they don't investigate further, I find their claim convincing.
    This reinforced to me that the quality of the training data is more important than the species it originates from.
    Strong performance in C. elegans and Arabadopsis also suggests that the human model is sufficient for diverse organisms.
    """
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
    # TODO: update rnasamba installation according to tests
    conda: "envs/rnasamba.yml"
    shell:'''
    rnasamba classify -p {output.fa} {output.tsv} {input.all_contigs} {input.model}
    '''

## TODO: predict sORFs from lncRNAs

################################################################################
## cleavage prediction
################################################################################

rule cleavage:
    """
    Defines a target rule for cleavage prediction so a user can run only cleavage prediction if they desire.
    snakemake cleavage --software-deployment-method conda -j 8 
    """
    input:
        "outputs/cleavage/nlpprecursor/nlpprecursor_ripp_predictions.tsv",
        "outputs/cleavage/deeppeptide/peptides.faa"

rule remove_stop_codon_asterisk_from_transdecoder_ORFs:
    input: orfs_amino_acids
    output: "outputs/cleavage/preprocessing/noasterisk.faa"
    shell:'''
    sed '/^[^>]/s/\*//g' {input} > {output}
    '''

# Ribosomally synthesized and post-translationally modified peptide prediction

rule download_nlpprecursor_models:
    output:
        tar = "inputs/models/nlpprecursor/nlpprecursor_models.tar.gz",
        model = "inputs/models/nlpprecursor/models/annotation/model.p"
    params: outdir = "inputs/models/nlpprecursor"
    shell:'''
    curl -JLo {output.tar} https://github.com/magarveylab/NLPPrecursor/releases/download/1.0/nlpprecursor_models.tar.gz
    tar xf {output.tar} -C {params.outdir} 
    '''

rule nlpprecursor:
    """
    The nlpprecursor tool is part of the [DeepRiPP approach](https://doi.org/10.1073/pnas.1901493116) that predicts ribosomally synthesized postranslationally modified peptides, a subclass of cleavage peptides.
    Unlike many tools in the RiPP prediction space, "NLPPrecursor identifies RiPPs independent of genomic context and neighboring biosynthetic genes."
    Note the paper suggests, "the precursor cleavage algorithm predicted N-terminal cleavage sites with 90% accuracy, when considering cleavage points ±5 amino acids from the true prediction site, a range within which all possible complete chemical structures can be elaborated in silico by combinatorial structure prediction."    
    From the DeepRiPP paper supplement: "Protein sequences of open reading frames are used as input.
    The output of the model consists of a classification of each ORF as either a precursor peptide (further subclassified according to RiPP family), or a non-precursor peptide. 
    A total of 14 classes are identified (n_class)."
    """
    input: 
        faa = "outputs/cleavage/preprocessing/noasterisk.faa",
        model = "inputs/models/nlpprecursor/models/annotation/model.p"
    output: "outputs/cleavage/nlpprecursor/nlpprecursor_ripp_predictions.tsv"
    params: modelsdir = "inputs/models/nlpprecursor/models/"
    conda: "envs/nlpprecursor.yml"
    shell:'''
    python scripts/run_nlpprecursor.py {params.modelsdir} {input.faa} {output}
    '''

# General Cleavage peptide prediction

rule clone_deeppeptide:
    output: "cloned_repositories/DeepPeptide/LICENSE"
    shell:'''
    cd cloned_repositories
    git clone https://github.com/fteufel/DeepPeptide.git
    '''

rule deeppeptide:
    input:
        src = "cloned_repositories/DeepPeptide/LICENSE",
        faa = "outputs/cleavage/preprocessing/noasterisk.faa"
    output: "outputs/cleavage/deeppeptide/peptide_predictions.json"
    conda: "envs/deeppeptide.yml"
    params: outdir = "outputs/cleavage/deeppeptide/"
    shell:'''
    cd cloned_repositories/DeepPeptide/predictor && python3 predict.py --fastafile ../../../{input.faa} --output_dir {params.outdir} --output_fmt json
    mv outputs/cleavage/deeppeptide/ ../../../outputs/cleavage/
    '''

rule extract_deeppeptide_sequences:
    """
    DeepPeptide outputs a json file of peptide predictions and locations, but does not output the sequences themselves.
    This step parses the JSON file and the protein FASTA from which peptides were predicted.
    It outputs the propeptide (full ORF, uncleaved) as well as the predicted peptide sequence (cleaved) in FASTA format.
    """
    input:  
        faa = "outputs/cleavage/preprocessing/noasterisk.faa",
        json = "outputs/cleavage/deeppeptide/peptide_predictions.json"
    output:
        propeptide = "outputs/cleavage/deeppeptide/propeptides.faa" ,
        peptide = "outputs/cleavage/deeppeptide/peptides.faa"
    conda: "envs/deeppeptide.yml"
    shell:'''
    python scripts/extract_deeppeptide_sequences.py {input.json} {input.faa} {output.propeptide} {output.peptide}
    '''
