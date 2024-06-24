from Bio import SeqIO
from Bio.Seq import Seq


#def verify_translation(nucleotide_seq, amino_acid_seq, to_stop):
#    """Verify that a nucleotide sequence translates correctly to its amino acid sequence."""
#    translated_seq = Seq(str(nucleotide_seq)).translate(to_stop=to_stop)
#    return str(translated_seq) == str(amino_acid_seq)

def verify_translation(nucleotide_seq, amino_acid_seq, to_stop, allow_wildcard_x = False):
    """
    Verify that a nucleotide sequence correctly translates to its corresponding amino acid sequence.
    The allow_wildcard_x argument permits any amino acid encoded by an "X" in the input sequence to
    be replaced with the corresponding translated nucleotide sequence. We include this due to
    Orfipy's conservative approach of not assigning an amino acid identity if the codon contains an
    "N". Conversely, Biopython is less restrictive and will assign an amino acid sequence if the
    codon can only translate into a single amino acid, regardless of the "N". For example, the
    leucine codon begins with CT. Therefore, even if the codon is CTN, we know that it encodes
    leucine. This approach requires the sequences to be the same length. Note that this function
    only verifies that two translations match; it does not edit the amino acid sequence. 
    """
    translated_seq = str(Seq(str(nucleotide_seq)).translate(to_stop=to_stop))
    amino_acid_seq = str(amino_acid_seq)

    if allow_wildcard_x:
        if len(translated_seq) != len(amino_acid_seq):
            return False
        amino_acid_seq = ''.join(
            translated_seq[i] if amino_acid == 'X' else amino_acid
            for i, amino_acid in enumerate(amino_acid_seq)
        )

    return translated_seq == amino_acid_seq


def read_fasta(fasta_file):
    """Read a FASTA file using BioPython and return a dictionary of sequences."""
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
    return sequences
