import os
from pathlib import Path

# This snakefile downloads and formats hidden markov model (HMM) profiles for non-ribosomal peptide
# synthetase (NRPS) annotation.
INPUT_DIR = Path("inputs/models/")
# Note input and output directories are the same.
# This snakefile sets up inputs for the main Snakefile.
# the output of this file is included in the repo ("inputs/models/nrps/nrps.hmm").
OUTPUT_DIR = Path("inputs/models/")
PFAM_NRPS_HMM_IDS = [
    "PF00501",
    "PF00975",
    "PF00668",
]


rule download_nrps_pfam_hmms:
    """
    This rule downloads PFAM HMM profiles for domains that are frequently annotated in NRPS enzymes.
    PFAM stands for protein families and is a database of protein families and domains.
    PFAM curates sequences and builds multiple sequence alignments and HMM profiles for annotations.
    We selected these sequences as examples from the literature frequently used to id NRPS enzymes.
    We provide examples citations below.
    An example of an NRPS with these domains is here: https://www.uniprot.org/uniprotkb/D2TVG0/entry
    PF00501 (AMP-binding enzyme; adenylation domain; 10.1016/j.ygeno.2022.110525)
    PF00975 (thioesterase domain https://www.rcsb.org/annotations/3FLB)
    PF00668 (condensation domain; c-domain superfamily members; 10.1073/pnas.2026017118)
    """
    output:
        hmm=INPUT_DIR / "nrps/pfam/{pfam_nrps_hmm_ids}.hmm",
    shell:
        """
        curl -JL https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{wildcards.nrps_pfam}?annotation=hmm | gunzip > {output.hmm}
        """


rule combine_and_convert_pfam_hmms:
    """
    This rule combines the PFAM domain HMM profiles.
    It then uses hmmconvert to format them for the version of HMMER used in this repo.
    """
    input:
        expand(rules.download_nrps_pfam_hmms.output.hmm, pfam_nrps_hmm_ids=PFAM_NRPS_HMM_IDS),
    output:
        hmm=OUTPUT_DIR / "nrps/pfam/pfam.hmm",
    conda:
        "envs/hmmer.yml"
    shell:
        """
        cat {input} | hmmconvert - > {output.hmm}
        """


rule download_antismash_hmms:
    """
    This rule downloads the NRPS and PKS domain HMM profiles used by the antismash annotation tool.
    It streams the file to hmmconvert to format the profiles for the version of HMMER in this repo.
    """
    output:
        hmm=OUTPUT_DIR / "nrps/antismash/nrpspksdomains.hmm",
    conda:
        "envs/hmmer.yml"
    shell:
        """
        curl -JL https://raw.githubusercontent.com/antismash/antismash/master/antismash/detection/nrps_pks_domains/data/nrpspksdomains.hmm | hmmconvert - > {output.hmm}
        """


rule download_nrps_motif_finder_hmms:
    """
    This rule downloads a C-domain HMM profiles from a recent PLOS study and combines.
    It formats them for the version of HMMER used in this repo.
    """
    output:
        hmm=OUTPUT_DIR / "nrps/nrps_motif_finder/nrps_c_domain.hmm",
    params:
        outdir=OUTPUT_DIR / "nrps/nrps_motif_finder/",
    conda:
        "envs/hmmer.yml"
    shell:
        """
        curl -JLo {params.outdir}/journal.pcbi.1011100.s056.zip https://doi.org/10.1371/journal.pcbi.1011100.s056
        unzip {params.outdir}/journal.pcbi.1011100.s056.zip -d {params.outdir}
        cat {params.outdir}/C\ domain\ subtype\ reference\ HMM\ files/*hmm | hmmconvert - > {output.hmm}
        # clean up intermediate files
        rm -rf {params.outdir}/journal.pcbi.1011100.s056.zip {params.outdir}/C\ domain\ subtype\ reference\ HMM\ files/
        """


rule combine_nrps_hmm_profiles:
    """
    This rule combines the three HMM profiles generated above into one file.
    This profile is used to identify NRPS domains in the main pipeline in this repo (Snakefile).
    These profiles are compatible because of the hmmconvert commands applied above.
    """
    input:
        rules.combine_and_convert_pfam_hmms.output.hmm,
        rules.download_antismash_hmms.output.hmm,
        rules.download_nrps_motif_finder_hmms.output.hmm,
    output:
        hmm=OUTPUT_DIR / "nrps/nrps.hmm",
    shell:
        """
        cat {input} > {output.hmm}
        """


rule all:
    default_target: True
    input:
        rules.combine_nrps_hmm_profiles.output.hmm,
