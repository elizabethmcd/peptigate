import argparse
import csv
import json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_fasta(fasta_file):
    """Read a FASTA file using BioPython and return a dictionary of sequences."""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences


def extract_peptide_sequences(
    data, fasta_file, proteins_output_file, peptides_output_file, predictions_output_file
):
    """
    Extract gene and peptide sequences based on the data dictionary and FASTA file,
    then write to separate files.

    Extracts protein and peptide sequences based on the provided data dictionary and a FASTA file.
    The data object is created from a JSON file output by deeppeptide.
    The protein sequences are extracted from the FASTA file using IDs in the data dictionary.
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
      Peptide sequences are also written in FASTA format, with headers indicating their source ID,
      start and end positions in the protein sequence, and that DeepPeptide was the source of the
      annotation.
    - predictions_output_file (str): Path to the output TSV file where predictions will be saved.

    Returns:
    None

    Raises:
    - FileNotFoundError: If the fasta_file does not exist or cannot be read.
    - KeyError: If the expected keys are not found in the data dictionary.

    Example usage:
    extract_peptide_sequences(
        data={'PREDICTIONS': {'>1': {'peptides': [{'start': 1, 'end': 9, 'type': 'Propeptide'}]}}},
        fasta_file='path/to/fasta_file.fasta',
        proteins_output_file='path/to/proteins_output.fasta',
        peptides_output_file='path/to/peptides_output.fasta',
        predictions_output_file='path/to/output.tsv')
    """
    sequences = read_fasta(fasta_file)

    protein_records = []
    peptide_records = []
    predictions = []

    for protein_key, protein_info in data["PREDICTIONS"].items():
        protein_id = protein_key.split()[0][1:]  # Extract the ID part
        peptides = protein_info.get("peptides")
        if peptides:  # Check if there are peptides
            protein_sequence = sequences.get(protein_id)
            if protein_sequence:  # If the protein sequence is found in the FASTA
                protein_records.append(
                    SeqRecord(Seq(protein_sequence), id=protein_id, description="")
                )

                for peptide in peptides:
                    start, end, peptide_class = peptide["start"], peptide["end"], peptide["type"]
                    peptide_metadata = {
                        "start": start,
                        "end": end,
                        "type": "cleavage",
                        "class": peptide_class,
                        "prediction_tool": "deeppeptide",
                    }
                    peptide_sequence = protein_sequence[start - 1 : end]  # Extract peptide sequence
                    peptide_id = f"{protein_id}_start{start}_end{end}"
                    description_fields = [
                        f"{key}:{value}" for key, value in peptide_metadata.items()
                    ]
                    peptide_records.append(
                        SeqRecord(
                            Seq(peptide_sequence),
                            id=peptide_id,
                            description=" ".join(description_fields),
                        )
                    )
                    predictions.append(
                        [peptide_id, start, end, "cleavage", peptide_class, "deeppeptide"]
                    )

    with open(proteins_output_file, "w") as proteins_out:
        SeqIO.write(protein_records, proteins_out, "fasta")

    with open(peptides_output_file, "w") as peptides_out:
        SeqIO.write(peptide_records, peptides_out, "fasta")

    with open(predictions_output_file, "w", newline="") as predictions_out:
        writer = csv.writer(predictions_out, delimiter="\t")
        writer.writerow(
            ["peptide_id", "start", "end", "peptide_type", "peptide_class", "prediction_tool"]
        )
        writer.writerows(predictions)


def main(
    json_file, fasta_file, proteins_output_file, peptides_output_file, predictions_output_file
):
    with open(json_file) as f:
        data = json.load(f)

    extract_peptide_sequences(
        data, fasta_file, proteins_output_file, peptides_output_file, predictions_output_file
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract peptide sequences from DeepPeptide JSON.")

    parser.add_argument("json_file", type=str, help="The JSON file output by DeepPeptide.")
    parser.add_argument("fasta_file", type=str, help="The protein FASTA file input to DeepPeptide.")
    parser.add_argument("proteins_output_file", type=str, help="The output file path for proteins.")
    parser.add_argument("peptides_output_file", type=str, help="The output file path for peptides.")
    parser.add_argument(
        "predictions_output_file", type=str, help="The output file path for predictions."
    )

    args = parser.parse_args()

    main(
        args.json_file,
        args.fasta_file,
        args.proteins_output_file,
        args.peptides_output_file,
        args.predictions_output_file,
    )
