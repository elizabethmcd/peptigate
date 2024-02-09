import json
import sys


def read_fasta(fasta_file):
    """Read a FASTA file and return a dictionary of sequences."""
    sequences = {}
    with open(fasta_file) as f:
        sequence_id = None
        sequence = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if sequence_id:
                    sequences[sequence_id] = sequence
                sequence_id = line[1:].split()[0]  # Assume ID is the first part after '>'
                sequence = ""
            else:
                sequence += line
        if sequence_id:
            sequences[sequence_id] = sequence
    return sequences


def extract_peptide_sequences(data, fasta_file, genes_output_file, peptides_output_file):
    """
    Extract gene and peptide sequences based on the data dictionary and FASTA file,
    then write to separate files.
    """
    sequences = read_fasta(fasta_file)

    with open(genes_output_file, "w") as genes_out, open(peptides_output_file, "w") as peptides_out:
        for transcript_key, transcript_info in data["PREDICTIONS"].items():
            transcript_id = transcript_key.split()[0][1:]  # Extract the ID part
            peptides = transcript_info.get("peptides", [])
            if peptides:  # Check if there are peptides
                gene_sequence = sequences.get(transcript_id)
                if gene_sequence:  # If the gene sequence is found in the FASTA
                    genes_out.write(f">{transcript_id}\n{gene_sequence}\n")
                    for peptide in peptides:
                        start, end = peptide["start"], peptide["end"]
                        peptide_sequence = gene_sequence[
                            start - 1 : end
                        ]  # Extract peptide sequence
                        peptides_out.write(
                            f">{transcript_id}_peptide_{start}_{end}\n{peptide_sequence}\n"
                        )


def main(json_file, fasta_file, genes_output_file, peptides_output_file):
    with open(json_file) as f:
        data = json.load(f)

    extract_peptide_sequences(data, fasta_file, genes_output_file, peptides_output_file)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python extract_deeppeptide_sequences.py <json> <fasta> <prot_out> <pep_out>")
        sys.exit(1)

    json_file = sys.argv[1]
    fasta_file = sys.argv[2]
    genes_output_file = sys.argv[3]
    peptides_output_file = sys.argv[4]

    main(json_file, fasta_file, genes_output_file, peptides_output_file)
