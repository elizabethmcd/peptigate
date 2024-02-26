import argparse
import csv


def read_deepsig_tsv(input_file):
    """
    Reads the TSV file output by deepsig and returns the data.

    Args:
        input_file (str): Path to the input deepsig TSV file.

    Returns:
        list of lists: Data from the TSV file.
    """
    with open(input_file, newline="") as file:
        reader = csv.reader(file, delimiter="\t")
        data = [row for row in reader]
    return data


def write_deepsig_tsv_with_header(output_file, data):
    """
    Adds headers to the deepsig data and writes to a TSV file.

    Args:
        output_file (str): Path to the output TSV file.
        data (list of lists): Data to be written to the file.
    """
    headers = [
        "peptide_id",
        "tool",
        "deepsig_feature",
        "deepsig_feature_start",
        "deepsig_feature_end",
        "deepsig_feature_score",
        "tmp1",
        "tmp2",
        "deepsig_description",
    ]

    with open(output_file, mode="w", newline="") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(headers)
        writer.writerows(data)


def main(input_file, output_file):
    """
    Reads and writes a deepsig TSV header, adding a header row.
    Deepsig does not include column names in its output.
    This script reads in the deepsig TSV, adds headers, and writes out the TSV.
    """
    data = read_deepsig_tsv(input_file)
    write_deepsig_tsv_with_header(output_file, data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add headers to a deepsig TSV file.")
    parser.add_argument("input_file", help="Path to the input deepsig TSV file.")
    parser.add_argument("output_file", help="Path to the output TSV file.")
    args = parser.parse_args()

    main(args.input_file, args.output_file)
