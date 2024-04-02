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

ORFS_AMINO_ACIDS = Path(config["orfs_amino_acids"])
ORFS_NUCLEOTIDES = Path(config["orfs_nucleotides"])
CONTIGS_SHORTER = Path(config["contigs_shorter_than_r2t_minimum_length"])
CONTIGS_LONGER = Path(config["contigs_longer_than_r2t_minimum_length"])
PLMUTILS_MODEL_DIR = Path(config["plmutils_model_dir"])

################################################################################
## Output file suffix conventions
################################################################################

"""
Note that we follow the conventions of the prokka tool for output file suffixes where possible.
* fna: nucleotide FASTA file of input contig sequences. These are transcripts from a transcriptome.
* faa: protein fasta file of translated CDS sequences.
    * faa: protein fasta file of translated CDS sequences for all predicted input proteins.
    * parent_faa: protein fasta file of translated CDS sequences for parent proteins of cleavage
                  peptides.
    * peptide_faa: protein fasta file of translated predicted peptide sequences
* ffn: nucleotide fasta file of CDS sequences.
    * ffn: nucleotide fasta file of all CDS sequences for all predicted input CDSs
    * parent_ffn: nucleotide fasta file of CDS sequences for parent CDSs of cleavage peptides
    * peptide_ffn: nucleotide fasta file of predicted peptide sequences
"""


################################################################################
## sORF prediction
################################################################################


rule combine_contigs:
    """
    By default we assume that files provided to this pipeline are reads2transcriptome outputs.
    Reads2transcriptome outputs two files that contain contigs.
    The first we refer to as contigs_shorter_than_r2t_minimum (or CONTIGS_SHORTER),
    which are transcripts output by assemblers that did not meet the r2t runs minimum contig length.
    The second we refer to as contigs_longer_than_r2t_minimum (or CONTIGS_LONGER) and are assembled
    transcripts that passed the isoform clustering and decontamination steps of r2t.
    We no longer need these pools of transcripts differentiated, so we combine them in this rule. 
    If your input transcriptome only has one file of assembled contigs, sequences only need to be
    supplied in contigs_longer_than_r2t_minimum_length.
    contigs_shorter_than_r2t_minimum_length can be an empty file.
    """
    input:
        contigs_shorter=CONTIGS_SHORTER,
        contigs_longer=CONTIGS_LONGER,
    output:
        all_contigs=OUTPUT_DIR / "sORF" / "contigs" / "all_input_contigs.fna",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        cat {input.contigs_shorter} {input.contigs_longer} > {output.all_contigs}
        """


rule get_coding_contig_names:
    """
    Extract amino acid contig names and remove everything after the first period, 
    which are isoform labels.
    This file will be used to select all contigs that DID NOT have a transdecoder-detected ORFs. 
    """
    input:
        ORFS_AMINO_ACIDS,
    output:
        names=OUTPUT_DIR / "sORF" / "contigs" / "orfs_amino_acid_names.txt",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit seq -n {input} | sed 's/[.].*$//' > {output}
        """


rule filter_contigs_to_no_predicted_ORF:
    """
    The peptigate pipeline takes as input transcripts (provided in two files) and predicted or
    annotated coding genes (provided in nucleotide and amino acid formats).
    This rule removes contigs with coding genes from the full transcript file.
    It assumes that the contig names are the same between to the two files (everything before the
    first period) and uses an inverted grep on the sequence names in the coding file to
    eliminate transcripts that contain protein-coding genes.
    It keeps all other transcripts, regardless of length, to investigate the presence of an sORF
    later in the pipeline.

    We expect that read2transcriptome will often be used to used to create input files for
    peptigate, or that transdecoder will be used to predict which transcripts contain protein-coding
    genes (the r2t pipeline also uses transdecoder to predict open reading frames (ORFs) from
    transcripts). By default, only ORFs that are longer than 100 amino acids are kept by
    transdecoder. The peptigate pipeline predicts peptides that are 100 amino acids or shorter.
    This rule eliminates transcripts that contained a transdecoder-predicted ORF.
    """
    input:
        fna=rules.combine_contigs.output.all_contigs,
        names=rules.get_coding_contig_names.output.names,
    output:
        fna=OUTPUT_DIR / "sORF" / "contigs" / "contigs_with_no_annotated_orf.fna",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit grep -v -f {input.names} {input.fna} -o {output.fna}
        """


rule download_uniref50_database:
    output:
        db=INPUT_DIR / "databases" / "uniref50.fasta.gz",
    shell:
        """
        curl -JLo {output} https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
        """


rule make_diamond_db_from_uniref50_database:
    input:
        db=rules.download_uniref50_database.output.db,
    output:
        db=INPUT_DIR / "databases" / "uniref50.dmnd",
    params:
        dbprefix=INPUT_DIR / "databases" / "uniref50",
    conda:
        "envs/diamond.yml"
    shell:
        """
        diamond makedb --in {input.db} -d {params.dbprefix}
        """


rule diamond_blastx_transcripts_against_uniref50_database:
    """
    This rule BLASTp's all contigs that did not have a predicted gene against uniref50 to check if
    the transcript looks like other coding sequences that are long (i.e., do not code for sORFs).
    We add this step after we observed that very fragmented transcriptome assemblies had over-
    prediction of sORF.
    When we BLASTp'd the sORFs against different databases, very few predicted sequences had hits.
    However, when we BLASTx'd the transcripts that gave rise to those sORFs, almost half had hits to
    non-sORF proteins (length > 100 amino acids).
    This supports the hypothesis that fragmented transcripts may be translated by `plmutils` in the wrong open
    reading frame which leads to over-prediction of sORFs.
    See https://github.com/Arcadia-Science/peptigate/issues/24 for more details.
    """
    input:
        db=rules.make_diamond_db_from_uniref50_database.output.db,
        fna=rules.filter_contigs_to_no_predicted_ORF.output.fna,
    output:
        tsv=OUTPUT_DIR / "sORF" / "filtering" / "0_blastx" / "matches.tsv",
    params:
        dbprefix=INPUT_DIR / "databases" / "uniref50",
    conda:
        "envs/diamond.yml"
    threads: 8
    shell:
        """
        diamond blastx -p {threads} -d {params.dbprefix} -q {input.fna} -o {output.tsv} \
            --header simple \
            --outfmt 6 qseqid sseqid full_sseq pident length qlen slen qcovhsp scovhsp mismatch gapopen qstart qend sstart send evalue bitscore
        """


rule filter_transcript_uniref50_hits:
    """
    Filter the uniref50 hits to those with a low evalue against long proteins (>100 amino acids).
    """
    input:
        tsv=rules.diamond_blastx_transcripts_against_uniref50_database.output.tsv,
    output:
        tsv=OUTPUT_DIR / "sORF" / "filtering" / "1_filtered_blastx" / "matches.tsv",
        txt=OUTPUT_DIR / "sORF" / "filtering" / "1_filtered_blastx" / "match_names.txt",
    conda:
        "envs/pandas.yml"
    shell:
        """
        python scripts/filter_transcript_uniref50_hits.py --input {input.tsv} \
            --output-blast {output.tsv} --output-names {output.txt}
        """


rule filter_no_predicted_ORF_contigs_to_no_uniref50_long_hits:
    """
    This rule filters the contigs that had no predicted ORF (defined by input files) to those that
    have not hits to long proteins in uniref50 (>100 amino acids).
    These transcripts will then be scanned for sORFs by plmutils.
    """
    input:
        fna=rules.filter_contigs_to_no_predicted_ORF.output.fna,
        names=rules.filter_transcript_uniref50_hits.output.txt,
    output:
        fna=OUTPUT_DIR
        / "sORF"
        / "filtering"
        / "contigs_with_no_predicted_orf_and_no_uniref50_blast_hit.fna",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit grep -v -f {input.names} {input.fna} -o {output.fna}
        """


rule plmutils_translate:
    """
    This rule takes input nucleotide transcripts, detects the longest open reading frame, and
    translates it into amino acid sequences.
    """
    input:
        rules.filter_no_predicted_ORF_contigs_to_no_uniref50_long_hits.output.fna,
    output:
        peptide_faa=OUTPUT_DIR / "sORF" / "plmutils" / "translated_contigs.faa",
    conda:
        "envs/plmutils.yml"
    shell:
        """
        plmutils translate --longest-only --output-filepath {output} {input}
        """


rule length_filter_plmutils_translate_output:
    input:
        rules.plmutils_translate.output.peptide_faa,
    output:
        peptide_faa=OUTPUT_DIR / "sORF" / "plmutils" / "translated_contigs_filtered.faa",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit seq --max-len 100 -o {output} {input}
        """


rule plmutils_embed:
    """
    This rule embeds amino acid sequences produced by plmutils translated into the embedding space
    of a protein large language model.
    For now, plmutils only supports ESM.
    We use the 8M model because it is fast and there is some evidence that it works well for other
    peptide prediction tasks (https://www.biorxiv.org/content/10.1101/2023.11.13.566825v2).
    The parameter --layer-ind -1 means to extract the embedding from the last layer of the model.
    """
    input:
        rules.length_filter_plmutils_translate_output.output.peptide_faa,
    output:
        npy=OUTPUT_DIR / "sORF" / "plmutils" / "embedded_contigs_filtered.npy",
    conda:
        "envs/plmutils.yml"
    shell:
        """
        plmutils embed --model-name esm2_t6_8M_UR50D \
            --layer-ind -1 \
            --output-filepath {output.npy} \
            {input}
        """


rule plmutils_predict:
    input:
        embeddings=rules.plmutils_embed.output.npy,
        peptide_faa=rules.length_filter_plmutils_translate_output.output.peptide_faa,
        model=PLMUTILS_MODEL_DIR,
    output:
        csv=OUTPUT_DIR / "sORF" / "plmutils" / "plmutils_predictions.csv",
    conda:
        "envs/plmutils.yml"
    shell:
        """
        plmutils predict --model-dirpath {input.model} \
            --embeddings-filepath {input.embeddings} \
            --fasta-filepath {input.peptide_faa} \
            --output-filepath {output.csv}
        """


rule extract_plmutils_predicted_peptides:
    input:
        csv=rules.plmutils_predict.output.csv,
        peptide_faa=rules.length_filter_plmutils_translate_output.output.peptide_faa,
    output:
        names=OUTPUT_DIR / "sORF" / "plmutils" / "plmutils_peptide_names.faa",
        peptide_faa=OUTPUT_DIR / "sORF" / "plmutils" / "plmutils_peptides.faa",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        grep "positive" {input.csv} | cut -d, -f1 > {output.names} 
        seqkit grep -f {output.names} {input.peptide_faa} -o {output.peptide_faa}
        """


rule extract_plmutils_predicted_peptides_as_nucleotides:
    input:
        fna=rules.filter_no_predicted_ORF_contigs_to_no_uniref50_long_hits.output.fna,
        peptide_faa=rules.extract_plmutils_predicted_peptides.output.peptide_faa,
    output:
        peptide_ffn=OUTPUT_DIR / "sORF" / "plmutils" / "plmutils_peptides.ffn",
    conda:
        "envs/biopython.yml"
    shell:
        """
        python  scripts/extract_plmutils_nucleotide_sequences.py \
            -n {input.fna} \
            -p {input.peptide_faa} \
            -o {output.peptide_ffn}
        """


################################################################################
## cleavage prediction
################################################################################


rule remove_stop_codon_asterisk_from_transdecoder_ORFs:
    input:
        ORFS_AMINO_ACIDS,
    output:
        faa=OUTPUT_DIR / "cleavage" / "preprocessing" / "noasterisk.faa",
    shell:
        """
        sed '/^[^>]/s/\*//g' {input} > {output}
        """


# Ribosomally synthesized and post-translationally modified peptide prediction


rule download_nlpprecursor_models:
    output:
        tar=INPUT_DIR / "models" / "nlpprecursor" / "nlpprecursor_models.tar.gz",
        model=INPUT_DIR / "models" / "nlpprecursor" / "models" / "annotation" / "model.p",
    params:
        outdir=INPUT_DIR / "models" / "nlpprecursor",
    shell:
        """
        curl -JLo {output.tar} https://github.com/magarveylab/NLPPrecursor/releases/download/1.0/nlpprecursor_models.tar.gz
        tar xf {output.tar} -C {params.outdir} 
        """


rule filter_protein_sequences_with_nonstandard_amino_acids:
    """
    The NLPPrecursor tool only uses the 20 standard amino acids.
    This script removes sequences with nonstandard amino acids in the sequence.
    """
    input:
        faa=rules.remove_stop_codon_asterisk_from_transdecoder_ORFs.output.faa,
    output:
        faa=OUTPUT_DIR / "cleavage" / "preprocessing" / "noasterisk_nononstandardaa.faa",
    conda:
        "envs/nlpprecursor.yml"
    shell:
        """
        python scripts/filter_protein_sequences_with_nonstandard_amino_acids.py --input {input} --output {output}
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
        faa=rules.filter_protein_sequences_with_nonstandard_amino_acids.output.faa,
        ffn=ORFS_NUCLEOTIDES,
        model=rules.download_nlpprecursor_models.output.model,
    output:
        parent_faa=OUTPUT_DIR / "cleavage" / "nlpprecursor" / "nlpprecursor_peptide_parents.faa",
        parent_ffn=OUTPUT_DIR / "cleavage" / "nlpprecursor" / "nlpprecursor_peptide_parents.ffn",
        peptide_faa=OUTPUT_DIR / "cleavage" / "nlpprecursor" / "nlpprecursor_peptides.faa",
        peptide_ffn=OUTPUT_DIR / "cleavage" / "nlpprecursor" / "nlpprecursor_peptides.ffn",
        tsv=OUTPUT_DIR / "cleavage" / "nlpprecursor" / "nlpprecursor_predictions.tsv",
    params:
        modelsdir=INPUT_DIR / "models" / "nlpprecursor" / "models/",
    conda:
        "envs/nlpprecursor.yml"
    shell:
        """
        python scripts/run_nlpprecursor.py \
            {params.modelsdir} \
            {input.faa} \
            {input.ffn} \
            {output.parent_faa} \
            {output.parent_ffn} \
            {output.peptide_faa} \
            {output.peptide_ffn} \
            {output.tsv}
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
        json=OUTPUT_DIR / "cleavage" / "deeppeptide" / "peptide_predictions.json",
    conda:
        "envs/deeppeptide.yml"
    params:
        outdir=OUTPUT_DIR / "cleavage" / "deeppeptide",
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
        json=rules.deeppeptide.output.json,
        faa=rules.remove_stop_codon_asterisk_from_transdecoder_ORFs.output.faa,
        ffn=ORFS_NUCLEOTIDES,
    output:
        parent_faa=OUTPUT_DIR / "cleavage" / "deeppeptide" / "deeppeptide_peptide_parents.faa",
        parent_ffn=OUTPUT_DIR / "cleavage" / "deeppeptide" / "deeppeptide_peptide_parents.ffn",
        peptide_faa=OUTPUT_DIR / "cleavage" / "deeppeptide" / "deeppeptide_peptides.faa",
        peptide_ffn=OUTPUT_DIR / "cleavage" / "deeppeptide" / "deeppeptide_peptides.ffn",
        tsv=OUTPUT_DIR / "cleavage" / "deeppeptide" / "deeppeptide_predictions.tsv",
    conda:
        "envs/biopython.yml"
    shell:
        """
        python scripts/extract_deeppeptide_sequences.py \
            {input.json} \
            {input.faa} \
            {input.ffn} \
            {output.parent_faa} \
            {output.parent_ffn} \
            {output.peptide_faa} \
            {output.peptide_ffn} \
            {output.tsv}
        """


################################################################################
## Combine peptide predictions
################################################################################


rule combine_peptide_faa_predictions:
    input:
        sorf=rules.extract_plmutils_predicted_peptides.output.peptide_faa,
        nlpprecursor=rules.nlpprecursor.output.peptide_faa,
        deeppeptide=rules.extract_deeppeptide_sequences.output.peptide_faa,
    output:
        peptide_faa=OUTPUT_DIR / "predictions" / "peptides.faa",
    shell:
        """
        cat {input} > {output.peptide_faa}
        """


rule convert_peptide_faa_to_tsv:
    input:
        peptide_faa=rules.combine_peptide_faa_predictions.output.peptide_faa,
    output:
        tsv=OUTPUT_DIR / "predictions" / "peptides_faa.tsv",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit fx2tab --only-id {input} -o {output}
        """


rule combine_peptide_ffn_predictions:
    input:
        sorf=rules.extract_plmutils_predicted_peptides_as_nucleotides.output.peptide_ffn,
        nlpprecursor=rules.nlpprecursor.output.peptide_ffn,
        deeppeptide=rules.extract_deeppeptide_sequences.output.peptide_ffn,
    output:
        peptide_ffn=OUTPUT_DIR / "predictions" / "peptides.ffn",
    shell:
        """
        cat {input} > {output.peptide_ffn}
        """


rule convert_peptide_ffn_to_tsv:
    input:
        peptide_ffn=rules.combine_peptide_ffn_predictions.output.peptide_ffn,
    output:
        tsv=OUTPUT_DIR / "predictions" / "peptides_ffn.tsv",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit fx2tab --only-id {input} -o {output}
        """


# TODO: decide whether to output ffn/faa of parent transcripts as a combined file.
#       I haven't needed it yet in subsequent analyses, so not outputting for now.
#       Note plmutils does not currently output parent transcript and would not have a parent
#       protein sequence.

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
        db=INPUT_DIR / "databases" / "peptipedia.fasta.gz",
    shell:
        """
        curl -JLo {output} https://osf.io/dzycu/download 
        """


rule make_diamond_db_from_peptipedia_database:
    input:
        db=rules.download_peptipedia_database.output.db,
    output:
        db=INPUT_DIR / "databases" / "peptipedia.dmnd",
    params:
        dbprefix=INPUT_DIR / "databases" / "peptipedia",
    conda:
        "envs/diamond.yml"
    shell:
        """
        diamond makedb --in {input.db} -d {params.dbprefix}
        """


rule diamond_blastp_peptide_predictions_against_peptipedia_database:
    input:
        db=rules.make_diamond_db_from_peptipedia_database.output.db,
        peptide_faa=rules.combine_peptide_faa_predictions.output.peptide_faa,
    output:
        tsv=OUTPUT_DIR / "annotation" / "peptipedia" / "blastp_matches.tsv",
    params:
        dbprefix=INPUT_DIR / "databases" / "peptipedia",
    conda:
        "envs/diamond.yml"
    shell:
        """
        diamond blastp -d {params.dbprefix} -q {input.peptide_faa} -o {output.tsv} --header simple \
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
        peptide_faa=rules.combine_peptide_faa_predictions.output.peptide_faa,
    output:
        tsv=OUTPUT_DIR / "annotation" / "deepsig" / "deepsig_predictions.tsv",
    conda:
        "envs/deepsig.yml"
    shell:
        """
        deepsig -f {input} -o {output}.tmp -k euk
        python scripts/add_header_to_deepsig_tsv.py {output}.tmp {output}
        """


rule characterize_peptides:
    input:
        peptide_faa=rules.combine_peptide_faa_predictions.output.peptide_faa,
    output:
        tsv=OUTPUT_DIR / "annotation" / "characteristics" / "peptide_characteristics.tsv",
    conda:
        "envs/peptides.yml"
    shell:
        """
        python scripts/characterize_peptides.py {input.peptide_faa} {output.tsv}
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


rule download_autopeptideml_models:
    output:
        # TER TODO: the authors of autopeptideml sent me these models.
        # They said they're working on uploading them.
        # Once they're available, I need to add a rule to download them and update the input here
        # to be the rules syntax
        # https://github.com/IBM/AutoPeptideML/issues/6#issuecomment-1989494559
        model=INPUT_DIR
        / "models"
        / "autopeptideml"
        / "HPO_NegSearch_HP"
        / "{autopeptideml_model_name}_1"
        / "apml_config.json",
    shell:
        """
        touch {output} # will become a curl command or something, depending on how models are packaged
        """


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
        peptide_faa=rules.combine_peptide_faa_predictions.output.peptide_faa,
        model=rules.download_autopeptideml_models.output.model,
    output:
        tsv=OUTPUT_DIR
        / "annotation"
        / "autopeptideml"
        / "autopeptideml_{autopeptideml_model_name}.tsv",
    params:
        modelsdir=INPUT_DIR / "models" / "autopeptideml" / "HPO_NegSearch_HP/",
    conda:
        "envs/autopeptideml.yml"
    shell:
        """
        python scripts/run_autopeptideml.py \
            --input_fasta {input.peptide_faa} \
            --model_folder {params.modelsdir}/{wildcards.autopeptideml_model_name}_1/ensemble \
            --model_name {wildcards.autopeptideml_model_name} \
            --output_tsv {output.tsv}
        """


rule combine_peptide_predictions:
    input:
        nlpprecursor=rules.nlpprecursor.output.tsv,
        deeppeptide=rules.extract_deeppeptide_sequences.output.tsv,
        plmutils=rules.plmutils_predict.output.csv,
        faa_tab=rules.convert_peptide_faa_to_tsv.output.tsv,
        ffn_tab=rules.convert_peptide_ffn_to_tsv.output.tsv,
    output:
        tsv=OUTPUT_DIR / "predictions" / "peptide_predictions.tsv",
    conda:
        "envs/tidyverse.yml"
    shell:
        """
        Rscript scripts/combine_peptide_predictions.R \
            --nlpprecursor_path {input.nlpprecursor} \
            --deeppeptide_path {input.deeppeptide} \
            --plmutils_path {input.plmutils} \
            --faa_tab_path {input.faa_tab} \
            --ffn_tab_path {input.ffn_tab} \
            --output_predictions_path {output.tsv}
        """


rule combine_peptide_annotations:
    input:
        autopeptideml=expand(
            rules.run_autopeptideml.output.tsv, autopeptideml_model_name=AUTOPEPTIDEML_MODEL_NAMES
        ),
        deepsig=rules.run_deepsig.output.tsv,
        peptipedia=rules.diamond_blastp_peptide_predictions_against_peptipedia_database.output.tsv,
        characteristics=rules.characterize_peptides.output.tsv,
    output:
        tsv=OUTPUT_DIR / "predictions" / "peptide_annotations.tsv",
    params:
        autopeptidemldir=OUTPUT_DIR / "annotation" / "autopeptideml/",
    conda:
        "envs/tidyverse.yml"
    shell:
        """
        Rscript scripts/combine_peptide_annotations.R \
            --autopeptideml_dir {params.autopeptidemldir} \
            --deepsig_path {input.deepsig} \
            --peptipedia_path {input.peptipedia} \
            --characteristics_path {input.characteristics} \
            --output_annotations_path {output.tsv}
        """


################################################################################
## Target rule all
################################################################################


rule all:
    default_target: True
    input:
        rules.combine_peptide_predictions.output.tsv,
        rules.combine_peptide_annotations.output.tsv,


rule predict_sORF:
    """
    Defines a target rule for sORF prediction so a user can run only sORF prediction.
    snakemake predict_sORF --software-deployment-method conda -j 8 
    """
    input:
        sorf=rules.extract_plmutils_predicted_peptides_as_nucleotides.output.peptide_ffn,


rule predict_cleavage:
    """
    Defines a target rule for cleavage prediction so a user can run only cleavage prediction.
    snakemake predict_cleavage --software-deployment-method conda -j 8 
    """
    input:
        rules.nlpprecursor.output.peptide_faa,
        rules.extract_deeppeptide_sequences.output.peptide_faa,
