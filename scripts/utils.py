from Bio import SeqIO
from Bio.Seq import Seq


def verify_translation(nucleotide_seq, amino_acid_seq, to_stop, allow_wildcard_x=False):
    """
    Verify that a nucleotide sequence correctly translates to its corresponding amino acid sequence.

    The `allow_wildcard_x` argument permits any amino acid encoded by an "X" in the input sequence
    to be replaced with the corresponding translated nucleotide sequence. This option is needed when
    `amino_acid_seq` was translated by a tool like Orfipy that does not translate codons that
    contain an "N". This is because Biopython (which is used here to translate `nucleotide_seq` into
    an amino acid sequence) is less restrictive and will translate codons containing "N"s if they
    can only translate into a single amino acid, regardless of the "N". For example, the leucine
    codons all begins with CT. Therefore, even if the codon is CTN, Biopython translates it to
    leucine.

    Notes:
    - This function requires the sequences to be the same length.
    - This function only verifies that two translations match; it does not edit the amino acid
      sequence.
    """
    translated_seq = str(Seq(str(nucleotide_seq)).translate(to_stop=to_stop))
    amino_acid_seq = str(amino_acid_seq)

    if len(translated_seq) != len(amino_acid_seq):
        return False

    if allow_wildcard_x:
        amino_acid_seq = "".join(
            translated_seq[i] if amino_acid == "X" else amino_acid
            for i, amino_acid in enumerate(amino_acid_seq)
        )

    return translated_seq == amino_acid_seq


def read_fasta(fasta_file):
    """Read a FASTA file using BioPython and return a dictionary of sequences."""
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
    return sequences
