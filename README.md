# peptigate 

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)
[![Snakemake](https://img.shields.io/badge/snakemake--green)](https://snakemake.readthedocs.io/en/stable/)

## Purpose

Peptigate investigates transcriptomes for peptides.
It predicts peptides (small open reading frames, cleavage peptides, and ribosomally translated and post-translationally modified) and then annotates them.

## Installation and Setup

This repository uses Snakemake to run the pipeline and conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/projects/miniconda/en/latest/). After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.

```{bash}
mamba env create -n peptigate --file envs/dev.yml
conda activate peptigate
```

Snakemake manages rule-specific environments via the `conda` directive and using environment files in the [envs/](./envs/) directory. Snakemake itself is installed in the main development conda environment as specified in the [dev.yml](./envs/dev.yml) file.

To start the pipeline, run:

```{bash}
snakemake --software-deployment-method conda -j 8
```

## Quick start

After setting up the Snakemake conda environment, clone the GitHub repository to your computer.
`cd` into the repository and use the `snakemake` command to run the pipeline on the demo data.

```
git clone https://github.com/Arcadia-Science/peptigate.git
cd peptigate
snakemake --software-deployment-method conda -j 1
```

## Input data

The [peptigate pipeline](./Snakefile) requires four input files:
* Open reading frames in amino acid format: Predicted open reading frames in amino acid format.
* Open reading frames in nucleotide format: Predicted open reading frames in nucleotide format. These open reading frames should have the same names as those in the amino acid file. Tools like [Transdecoder](https://github.com/TransDecoder/TransDecoder) provide these files in the correct format.
* Long transcripts/contigs: A transcriptome assembly FASTA file in nucleotide format containing transcripts or contigs.
* Short transcripts/contigs: contigs that are shorter than X nucleotides. Some transcriptome assemblers discard very short contigs and do not include them in the final assembly. However, some provide them as an intermediate output file. These contigs may contain sORFs and so are included as an input to the peptigate pipeline. If you do not have a file that contains very short contigs, provide a path to an empty file. Very short contigs can also be provided as part of the previous file. If that is the case, provide a path to an empty file for this input file.

The pipeline also requires pointers to three directories:
* Input directory path: This folder is used by the pipeline to store databases and models.
* Output directory path: This folder will be created by the snakemake pipeline and will be used to store the output files from the pipeline.
* Path to directory with plm-utils model: The path to the directory that stores the model for the plm-utils tool. A plm-utils model is provided in this repository.

These inputs are provided to the peptigate pipeline by a config file.
[`config.yml`](./config.yml) is an example config file, while the [`demo`](./demo) directory contains example input files.

We also included a [workflow](./protein_as_input.snakefile) that can take a single file of protein sequences as input.
[`config_protein.yml`](./config_protein.yml) is an example config file for this workflow.

Many of the steps in the workflow require models or databases.
These data are either included in the [`inputs`](./inputs) folder in this repository or are downloaded by the snakemake pipeline itself.

## Overview

Peptigate is a workflow that predicts bioactive peptides from transcriptome assemblies and annotates those predictions.
This Snakemake-based pipeline integrates tools to identify small open reading frames (sORFs) and cleavage peptides.
Each peptide prediction is then annotated to provide clues as to the bioactivity or function of the peptide.

### Description of how the tool works

Peptigate is broken into three sections.

#### Small open reading frame (sORF) prediction

**Background** Small open reading frames are peptides that are ["born small"](https://doi.org/10.1111/febs.15769) --
they are less than 100 amino acids but are otherwise like longer proteins in that they are synthesized via DNA transcription and ribosomal translation.
Many tools that predict open reading frames (ORFs) in transcripts have decreased accuracy at shorter lengths and by default do not output predictions shorter than 100 amino acids.
Yet, some sORFs produce functional proteins, leading to a systematic underappreciation in the detection and biological role of these proteins.
Techniques like ribosomal profiling and peptidomics mass spectrometry have highlighted that ubiquity of these proteins as well as some of their biological roles [https://doi.org/10.3389/fgene.2021.796060].
Many sORFs are located upstream or downstream of canonical long ORFs and play a regulatory role by influencing translation [https://doi.org/10.3389/fgene.2021.796060].
Other sORFs are encode functional peptides [https://doi.org/10.3389/fgene.2021.796060; https://doi.org/10.1021/pr401280w].

**What the pipeline does** The peptigate pipeline targets stand-alone sORFs with the goal of identifying functional peptides.
Peptigate begins sORF prediction by removing transcripts that had predicted ORFs.
It then uses the [`plm-utils`](https://github.com/Arcadia-Science/2024-plm-utils) tool to predict open reading frames from the rest of the transcripts.
plm-utils uses the start codons [TTG, CTG, ATG, GTG, ACG](https://doi.org/10.3389/fgene.2021.796060) and traditional stop codons for ORF prediction as sORFs frequently use non-canonical start sites.
If a predicted ORF is shorter than 301 nucleotides, it then predicts whether the ORF is coding or not using a model trained on [ESM embeddings](https://github.com/facebookresearch/esm).
We think these embeddings capture information about the secondary structure of proteins; large protein language models have [previously been shown](https://doi.org/10.1016/j.str.2022.11.012) to have high accuracy on structural predictions for some peptides. 
Peptigate returns the peptide sequences of predicted peptides in amino acid and nucleotide format.

**Observations from running the pipeline**: Because many valid sORFs are co-encoded on transcripts with longer ORFs, peptigate will over-predict functional sORFs when fragmented transcripts are supplied to the pipeline.
The peptigate will detect sORFs encoded in the 5' or 3' UTR of the longer canonical ORF if the canonical ORF was not annotated because it is part of a fragmented contig.
We implemented a filtering step to try and catch some of these cases but this is a tricky problem for which we don't have a great solution.
If this is your situation, two tricks that might prove useful to filter predictions down to the most likely functional peptides are to:
1. Search for sORF predictions with both chain and signal peptides annotated. If a peptide has both, it may be more likely to be targetted to a specific location, potentially indicating that is is more likely to be functional.
2. Filter to sORF predictions that have hits against the peptipedia database. This will limit to peptide predictions that are homologous to peptides that have been discovered before. 

#### Cleavage peptide prediction

#### Predicted peptide annotation

### Description of the folder structure and files

#### Folders and files in this repository

* [demo](./demo): Small input files for a demo run of the pipeline. If you're unsure what format your input files need to be in, cross-check against files in the `demo` data folder.
* [envs](./envs): This repository uses conda to manage software installations and versions. All software required for peptigate use and development is recorded in this folder.
* [inputs](./inputs/models): Required inputs for the peptigate pipeline. Currently includes the plmutils (sORF prediction) model as well as the pointers to the data used to train this model.
* [scripts](./scripts): Python and R scripts used by the Snakefiles in this repository.
* [`LICENSE`](./LICENSE): License specifying the re-use terms for the code in this repository.
* [`README.md`](./README.md): File outlining the contents of this repository and how to use the peptigate tool.
* [`Snakefile`](./Snakefile): The snakemake workflow file that orchestrates the full peptigate pipeline. 
* [`config.yml`](./config.yml): A demo & template config file for the main `Snakefile`.
* [`protein_as_input.snakefile`](./protein_as_input.snakefile): A simplified version of the peptigate pipeline that only requires a FASTA file of proteins (amino acid format).
* [`config_protein.yml`](./config_protein.yml): A demo & template config file for the `protein_as_input.snakefile`.
* [`curate_datasets_and_build_models.snakefile`](./curate_datasets_and_build_models.snakefile): workflow recording how we generated the [plm-utils sORF prediction model](./inputs/models/plmutils/). Because we provide this model in the inputs folder of this repository, we do not anticipate that most users will be interested in running this workflow. 
* [.github](./.github), [.vscode](./.vscode), [Makefile](./Makefile), [pyproject.toml](./Makefile): Control the developer behavior of the repository. See the [template repository](https://github.com/Arcadia-Science/snakemake-template) for a description of how these files work.

#### Folders and files output by the workflow

All predicted peptide sequences and annotation information are reported in the `predictions` folder.
Other folders record intermediate files needed to make the final prediction files.
See below for a description of each folder.
 
* annotation
    * autopeptideml: intermediate files recording bioacitivty predictions made by AutoPeptideML models.
    * characteristics: intermediate files recording chemical characteristics of predicted peptides calculated by the Python package peptides.
    * deepsig: intermediate files recording signal peptide predictions for the predicted peptides calculated by DeepSig.
    * peptipedia: intermediate files record BLASTp matches against the Peptipedia database using diamond.
* cleavage
    * deeppeptide: intermediate files and cleavage peptide predictions made by DeepPeptide
    * nlpprecursor: intermediate files and cleavage peptide predictions made by NLPPrecursor.
    * preprocessing: intermediate files formatted in preparation for cleavage peptide prediction.
* predictions: combined peptide predictions and annotations.
    * `peptide_annotations.tsv`: annotations for peptide predictions. Note a single `peptide_id` may have multiple rows, usually because DeepSig predicts both a chain and a signal peptide in a peptide sequence. 
    * `peptide_predictions.tsv`: peptide predictions, including protein and nucleotide sequences and prediction tool.
    * `peptides.faa`: peptide predictions in amino acid FASTA format.
    * `peptides.fna`: peptide predictions in nucleotide FASTA format.
    * `peptides_faa.tsv`: peptide predictions in amino acid format as a TSV file.
    * `peptides_fna.tsv`: peptide predictions in nucleotide format as a TSV file.
* sORF: intermediate files related to sORF prediction from non-coding nucleotide contigs supplied as input to the pipeline.
    * contigs: contains processed contiguous sequences that are filtered before sORF prediction.
    * filtering: a BLAST filter designed to remove fragmented contiguous sequences that contain fragmented coding sequences (canonical genes longer than 100 amino acids).
    * plmutils: intermediate files and sORF prediction via the tool plm-utils.

### Compute Specifications

We ran the pipeline on an AWS EC2 instance type `g4dn.2xlarge` running AMI Deep Learning Base OSS Nvidia Driver GPU AMI (Ubuntu 20.04) 20240122 (AMI ID ami-07eb000b3340966b0).
The tools plm-utils, DeepPeptide, NLPPrecursor, and AutoPeptideML can use GPUs so compute times will be substantially faster (hours vs. days) on a GPU than a CPU.
We did not test the pipeline on a CPU so there is a chance it will only work on a GPU.

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide-credit-for-contributions.md).
