"""
This script filters FASTA files to include only sequences containing specific amino acids.
It was inspired by the 20 amino acid library implemented in the NLPPrecursor models.

Usage:
    python filter_protein_sequences_with_nonstandard_amino_acids.py \
        --input INPUT_FASTA_FILE \
        --output OUTPUT_FASTA_FILE

Arguments:
    --input: Path to the input FASTA file containing amino acid sequences.
    --output: Path to the output FASTA file to save filtered sequences.

The script filters out any sequences that contain amino acids outside of the following allowed set:
M, A, N, F, D, L, V, Q, K, S, G, Y, T, C, R, H, I, P, E, W, X.
"""

import argparse

from Bio import SeqIO


def is_allowed_sequence(seq, allowed_amino_acids):
    """
    Check if all amino acids in a sequence are within the allowed set.

    Parameters:
    - seq (str): The amino acid sequence to check.
    - allowed_amino_acids (set): A set of allowed amino acid single-letter codes.

    Returns:
    - bool: True if the sequence contains only allowed amino acids, False otherwise.
    """
    return all(amino_acid in allowed_amino_acids for amino_acid in seq)


def filter_fasta(input_file, output_file, allowed_amino_acids):
    """
    Filter sequences in a FASTA file based on allowed amino acids.

    Parameters:
    - input_file (str): Path to the input FASTA file.
    - output_file (str): Path to the output FASTA file where filtered sequences are saved.
    - allowed_amino_acids (set): A set of allowed amino acid single-letter codes.
    """
    with open(output_file, "w") as filtered:
        for record in SeqIO.parse(input_file, "fasta"):
            if is_allowed_sequence(str(record.seq), allowed_amino_acids):
                SeqIO.write(record, filtered, "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter FASTA files to include only sequences with specific amino acids."
    )
    parser.add_argument("--input", required=True, help="Path to the input FASTA file.")
    parser.add_argument("--output", required=True, help="Path to the output FASTA file.")

    args = parser.parse_args()

    # This set represents the 20 standard amino acids.
    allowed_amino_acids = {
        "M",
        "A",
        "N",
        "F",
        "D",
        "L",
        "V",
        "Q",
        "K",
        "S",
        "G",
        "Y",
        "T",
        "C",
        "R",
        "H",
        "I",
        "P",
        "E",
        "W",
        "X",
    }

    filter_fasta(args.input, args.output, allowed_amino_acids)
