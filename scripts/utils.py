from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def verify_translation(nucleotide_seq, amino_acid_seq, to_stop):
    """Verify that a nucleotide sequence translates correctly to its amino acid sequence."""
    translated_seq = Seq(nucleotide_seq).translate(to_stop=to_stop)
    return str(translated_seq) == amino_acid_seq

def read_fasta(fasta_file):
    """Read a FASTA file using BioPython and return a dictionary of sequences."""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences
