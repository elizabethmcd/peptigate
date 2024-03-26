import argparse
import csv
import json

import utils
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def extract_peptide_sequences(
    data,
    protein_fasta_file,
    nucleotide_fasta_file,
    proteins_output_file,
    nucleotides_output_file,
    protein_peptides_output_file,
    nucleotide_peptides_output_file,
    predictions_output_file,
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
    - protein_fasta_file (str): The path to a FASTA file containing protein sequences.
      This should be the same file used to make the DeepPeptide predictions.
    - nucleotide_fasta_file (str): The path to a FASTA file containing nucleotide sequences.
      This should include sequences for the same genes as used to make the DeepPeptide predictions.
    - proteins_output_file (str): The path to the output file where protein sequences that gave rise
      to predicted peptides will be saved. Each sequence is written in amino acid FASTA format with
      its ID as the header.
    - nucleotides_output_file (str): The path to the output file where gene sequences that gave rise
      to predicted peptides will be saved. Each sequence is written in FASTA format with its ID as
      the header.
    - protein_peptides_output_file (str): The path to the output file where peptide sequences will
      be saved. Peptide sequences written to this file are saved in amino acid FASTA format, with
      headers indicating their source ID, start and end positions in the protein sequence, and that
      DeepPeptide was the source of the annotation.
    - nucleotide_peptides_output_file (str): The path to the output file where peptide sequences
      will be saved. Peptide sequences written to this file are saved in nucleotide FASTA format,
      with headers indicating their source ID, start and end positions in the protein sequence, and
      that DeepPeptide was the source of the annotation.
    - predictions_output_file (str): Path to the output TSV file where predictions will be saved.

    Returns:
    None

    Raises:
    - FileNotFoundError: If the fasta_file does not exist or cannot be read.
    - KeyError: If the expected keys are not found in the data dictionary.

    Example usage:
    extract_peptide_sequences(
        data={'PREDICTIONS': {'>1': {'peptides': [{'start': 1, 'end': 9, 'type': 'Propeptide'}]}}},
        protein_fasta_file='path/to/protein_fasta_file.faa',
        nucleotide_fasta_file='path/to/nucleotides_fasta_file.fna',
        proteins_output_file='path/to/proteins_output.faa',
        nucleotides_output_file='path/to/nucleotides_output.fna',
        protein_peptides_output_file='path/to/peptides_output.faa',
        nucleotide_peptides_output_file='path/to/peptides_output.fna',
        predictions_output_file='path/to/output.tsv')
    """
    protein_sequences = utils.read_fasta(protein_fasta_file)
    nucleotide_sequences = utils.read_fasta(nucleotide_fasta_file)

    protein_records = []
    nucleotide_records = []
    protein_peptide_records = []
    nucleotide_peptide_records = []
    predictions = []

    for protein_key, protein_info in data["PREDICTIONS"].items():
        protein_id = protein_key.split()[0][1:]

        peptides = protein_info.get("peptides")
        if peptides:
            protein_sequence = protein_sequences.get(protein_id)
            nucleotide_sequence = nucleotide_sequences.get(protein_id)
            if protein_sequence:
                protein_records.append(
                    SeqRecord(Seq(protein_sequence), id=protein_id, description="")
                )
                nucleotide_records.append(
                    SeqRecord(Seq(nucleotide_sequence), id=protein_id, description="")
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
                    protein_peptide_sequence = protein_sequence[start - 1 : end]
                    nucleotide_peptide_sequence = nucleotide_sequence[(start - 1) * 3 : end * 3]
                    if utils.verify_translation(
                        nucleotide_peptide_sequence, protein_peptide_sequence, to_stop=False
                    ):
                        peptide_id = f"{protein_id}_start{start}_end{end}"
                        description_fields = [
                            f"{key}:{value}" for key, value in peptide_metadata.items()
                        ]
                        description = " ".join(description_fields)
                        protein_peptide_records.append(
                            SeqRecord(
                                Seq(protein_peptide_sequence),
                                id=peptide_id,
                                description=description,
                            )
                        )
                        nucleotide_peptide_records.append(
                            SeqRecord(
                                Seq(nucleotide_peptide_sequence),
                                id=peptide_id,
                                description=description,
                            )
                        )
                        predictions.append(
                            [peptide_id, start, end, "cleavage", peptide_class, "deeppeptide"]
                        )

    SeqIO.write(protein_records, proteins_output_file, "fasta")
    SeqIO.write(nucleotide_records, nucleotides_output_file, "fasta")
    SeqIO.write(protein_peptide_records, protein_peptides_output_file, "fasta")
    SeqIO.write(nucleotide_peptide_records, nucleotide_peptides_output_file, "fasta")

    with open(predictions_output_file, "w", newline="") as predictions_out:
        writer = csv.writer(predictions_out, delimiter="\t")
        writer.writerow(
            ["peptide_id", "start", "end", "peptide_type", "peptide_class", "prediction_tool"]
        )
        writer.writerows(predictions)


def main(
    json_file,
    protein_fasta_file,
    nucleotide_fasta_file,
    proteins_output_file,
    nucleotides_output_file,
    protein_peptides_output_file,
    nucleotide_peptides_output_file,
    predictions_output_file,
):
    with open(json_file) as f:
        data = json.load(f)

    extract_peptide_sequences(
        data,
        protein_fasta_file,
        nucleotide_fasta_file,
        proteins_output_file,
        nucleotides_output_file,
        protein_peptides_output_file,
        nucleotide_peptides_output_file,
        predictions_output_file,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract peptide sequences from DeepPeptide JSON.")

    parser.add_argument("json_file", type=str, help="The JSON file output by DeepPeptide.")
    parser.add_argument(
        "protein_fasta_file", type=str, help="The protein FASTA file input to DeepPeptide."
    )
    parser.add_argument(
        "nucleotide_fasta_file", type=str, help="The nucleotide FASTA file for the genes."
    )
    parser.add_argument("proteins_output_file", type=str, help="The output file path for proteins.")
    parser.add_argument(
        "nucleotides_output_file", type=str, help="The output file path for nucleotide sequences."
    )
    parser.add_argument(
        "protein_peptides_output_file",
        type=str,
        help="The output file path for peptides in amino acid format.",
    )
    parser.add_argument(
        "nucleotide_peptides_output_file",
        type=str,
        help="The output file path for peptides in nucleotide format.",
    )
    parser.add_argument(
        "predictions_output_file", type=str, help="The output file path for predictions."
    )

    args = parser.parse_args()

    main(
        args.json_file,
        args.protein_fasta_file,
        args.nucleotide_fasta_file,
        args.proteins_output_file,
        args.nucleotides_output_file,
        args.protein_peptides_output_file,
        args.nucleotide_peptides_output_file,
        args.predictions_output_file,
    )
