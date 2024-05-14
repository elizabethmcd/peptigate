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

## Data

The peptigate pipeline requires four input files:

Only have predicted proteins, perhaps from a genome or from some other source?


## Overview

### Description of how the tool works

### Description of the folder structure and files

#### Folders and files in this repository

* [demo](./demo):
* [envs](./envs):
* [inputs](./inputs/models):
* [scripts](./scripts):
* [LICENSE](./LICENSE): License specifying the re-use terms for the code in this repository.
* [README.md](./README.md):
* [Snakefile](./Snakefile):
* [config.yml](./config.yml):
* [config_protein.yml](./config_protein.yml)
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
The tools plm-utils, DeepPeptide, NLPPrecursor, and AutoPeptideML all rely on a GPU so compute times will be substantially faster (hours vs. days) on a GPU than a CPU.
We did not test the pipeline on a CPU so there is a chance it will only work on a GPU.

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide-credit-for-contributions.md).
