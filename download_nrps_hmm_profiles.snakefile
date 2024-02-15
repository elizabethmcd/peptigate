import os
from pathlib import Path

INPUT_DIR = Path("inputs/models/")
OUTPUT_DIR = Path("inputs/models/")
NRPS_PFAMS = ['PF00501', 'PF06339', 'PF00975', 'PF00668',]

rule download_nrps_pfam_hmms:
    output: hmm = INPUT_DIR / "nrps/pfam/{nrps_pfam}.hmm"
    shell:
        """
        curl -JL https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{wildcards.nrps_pfam}?annotation=hmm | gunzip > {output.hmm}
        """

rule combine_and_convert_pfam_hmms:
    input: expand(rules.download_nrps_pfam_hmms.output.hmm, nrps_pfam = NRPS_PFAMS)
    output: hmm=OUTPUT_DIR / "nrps/pfam/pfam.hmm"
    conda: "envs/hmmer.yml"
    shell:
        """
        cat {input} | hmmconvert - > {output.hmm}
        """

rule download_antismash_hmms:
    output: hmm = OUTPUT_DIR / "nrps/antismash/nrpspksdomains.hmm"
    conda: "envs/hmmer.yml"
    shell:
        """
        curl -JL https://raw.githubusercontent.com/antismash/antismash/master/antismash/detection/nrps_pks_domains/data/nrpspksdomains.hmm | hmmconvert - > {output.hmm}
        """

rule download_nrps_motif_finder_hmms:
    output:
        archive = INPUT_DIR / "nrps/nrps_motif_finder/journal.pcbi.1011100.s056.zip",
        hmm = OUTPUT_DIR / "nrps/nrps_motif_finder/nrps_c_domain.hmm"
    params: outdir = OUTPUT_DIR / "nrps/nrps_motif_finder/"
    conda: "envs/hmmer.yml"
    shell:
        """
        curl -JLo {output.archive} https://doi.org/10.1371/journal.pcbi.1011100.s056
        unzip {output.archive} -d {params.outdir}
        cat {params.outdir}/C\ domain\ subtype\ reference\ HMM\ files/*hmm | hmmconvert - > {output.hmm}
        """

rule combine_nrps_hmm_profiles:
    input:
        rules.combine_and_convert_pfam_hmms.output.hmm,
        rules.download_antismash_hmms.output.hmm,
        rules.download_nrps_motif_finder_hmms.output.hmm
    output: hmm=OUTPUT_DIR / "nrps/nrps.hmm"
    shell:
        """
        cat {input} > {output.hmm}
        """

rule all:
    default_target: True
    input:
       rules.combine_nrps_hmm_profiles.output.hmm
