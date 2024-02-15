import argparse
import json
import sys
from Bio import SeqIO

def read_fasta(fasta_file):
    """Read a FASTA file using BioPython and return a dictionary of sequences."""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences


def extract_peptide_sequences(data, fasta_file, proteins_output_file, peptides_output_file):
    """
    Extract gene and peptide sequences based on the data dictionary and FASTA file,
    then write to separate files.

    Extracts protein and peptide sequences based on the provided data dictionary and a FASTA file.
    The data object is created from a JSON file output by deeppeptide.
    The protein sequences are extracted from the FASTA file using the IDs found in the data dictionary.
    Peptide sequences are then extracted from these protein sequences based on start and end 
    positions specified for each peptide within the data dictionary. 
    The extracted protein and peptide sequences are written to separate output files.

    Parameters:
    - data (dict): A dictionary containing prediction data, where each key is a protein ID and 
      the associated value is another dictionary with details including peptides' start and end 
      positions.
    - fasta_file (str): The path to a FASTA file containing protein sequences.
      This should be the same file used to make the DeepPeptide predictions.
    - proteins_output_file (str): The path to the output file where protein sequences will be saved.
      Each sequence is written in FASTA format with its ID as the header.
    - peptides_output_file (str): The path to the output file where peptide sequences will be saved.
      Peptide sequences are also written in FASTA format,
      with headers indicating their source transcript ID and their start and end positions within
      the protein sequence.

    Returns:
    None

    Raises:
    - FileNotFoundError: If the fasta_file does not exist or cannot be read.
    - KeyError: If the expected keys are not found in the data dictionary.

    Example usage:
    extract_peptide_sequences(data={'PREDICTIONS': {'>1': {'peptides': [{'start': 1, 'end': 9}]}}},
                              fasta_file='path/to/fasta_file.fasta',
                              proteins_output_file='path/to/proteins_output.fasta',
                              peptides_output_file='path/to/peptides_output.fasta')
    """
    sequences = read_fasta(fasta_file)

    with open(proteins_output_file, "w") as proteins_out, open(peptides_output_file, "w") as peptides_out:
        for protein_key, protein_info in data["PREDICTIONS"].items():
            protein_id = protein_key.split()[0][1:]  # Extract the ID part
            peptides = protein_info.get("peptides", [])
            if peptides:  # Check if there are peptides
                protein_sequence = sequences.get(protein_id)
                if protein_sequence:  # If the protein sequence is found in the FASTA
                    proteins_out.write(f">{protein_id}\n{protein_sequence}\n")
                    for peptide in peptides:
                        start, end = peptide["start"], peptide["end"]
                        peptide_sequence = protein_sequence[
                            start - 1 : end
                        ]  # Extract peptide sequence
                        peptides_out.write(
                            f">{protein_id}_peptide_{start}_{end}\n{peptide_sequence}\n"
                        )


def main(json_file, fasta_file, proteins_output_file, peptides_output_file):
    with open(json_file) as f:
        data = json.load(f)

    extract_peptide_sequences(data, fasta_file, proteins_output_file, peptides_output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract peptide sequences from DeepPeptide JSON.')

    # Add the arguments
    parser.add_argument('json_file', type=str, help='The JSON file output by DeepPeptide.')
    parser.add_argument('fasta_file', type=str, help='The protein FASTA file input to DeepPeptide.')
    parser.add_argument('proteins_output_file', type=str, help='The output file path for proteins.')
    parser.add_argument('peptides_output_file', type=str, help='The output file path for peptides.')

    # Execute the parse_args() method
    args = parser.parse_args()

    main(args.json_file, args.fasta_file, args.proteins_output_file, args.peptides_output_file)
