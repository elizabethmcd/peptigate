import argparse

import pandas as pd


def filter_blast_results(input_file, output_blast, output_names):
    """
    Filters DIAMOND BLAST results to identify and retain hits that are likely to correspond to long
    protein-coding sequences, based on predefined criteria of e-value and sequence length (slen).

    This function filters out BLAST hits with an e-value less than 1e-10, a subject (hit)
    sequence length (slen) of less than 100 amino acids, and an alignment length of greater than 100
    amino acids. These thresholds are chosen to select transcripts likely hitting long protein
    coding sequences in databases like UniRef50, thus suggesting that these transcripts are less
    likely to encode short open reading frame (sORF) peptides because they are likely fragments of
    true long protein-coding genes.

    Parameters:
    - input_file (str): Path to the input file containing DIAMOND BLAST results.
    - output_blast (str): Path to the output file where the filtered BLAST results will be saved.
    - output_names (str): Path to the output file where the qseqid of retained hits will be saved.

    Outputs are saved as tab-separated values (TSV) files.
    The first output file contains all columns of the filtered BLAST results,
    and the second output file contains only the qseqid of the retained hits.
    """
    blast_hits = pd.read_csv(input_file, sep="\t")
    filtered_blast_hits = blast_hits[
        (blast_hits["evalue"] < 1e-10) & (blast_hits["slen"] > 100) & (blast_hits["length"] > 100)
    ]
    filtered_sorted_blast_hits = filtered_blast_hits.sort_values(
        by=["qseqid", "bitscore", "evalue"], ascending=[True, False, True]
    )
    final_blast_hits = filtered_sorted_blast_hits.drop_duplicates(subset="qseqid", keep="first")

    # Save the filtered BLAST results to a new file.
    final_blast_hits.to_csv(output_blast, sep="\t", index=False)

    # Save a file with just the qseqid of retained hits.
    final_blast_hits["qseqid"].to_csv(output_names, index=False, header=False)


def main():
    parser = argparse.ArgumentParser(description="Filter DIAMOND BLAST results.")
    parser.add_argument("--input", help="Input file path for DIAMOND BLAST results.", required=True)
    parser.add_argument(
        "--output-blast", help="Output file path for filtered results.", required=True
    )
    parser.add_argument(
        "--output-names",
        help="Output file path for names of sequences in filtered results.",
        required=True,
    )

    args = parser.parse_args()

    filter_blast_results(args.input, args.output_blast, args.output_names)


if __name__ == "__main__":
    main()
