import csv
import sys
import time
from pathlib import Path
import argparse

import nlpprecursor
from Bio import SeqIO
from Bio.Seq import Seq
from nlpprecursor.annotation.data import DatasetGenerator as ADG
from nlpprecursor.classification.data import DatasetGenerator as CDG

# This allows for backwards compatibility of the pickled models.
sys.modules["protai"] = nlpprecursor


def robust_predict(predict_function, *args, max_attempts=2, sleep_time=1):
    """
    Attempts to call the predict function up to a maximum number of attempts.
    TODO: debug this problem if this pipeline becomes widely used
    On first attempt to call this function on a GPU,
    the function produces a RuntimeError: cuDNN error: CUDNN_STATUS_EXECUTION_FAILED.
    Second execution succeeds.
    Execution fails on the @classmethod predict(), on the line `predictions = model(tokens)[0]`.
    Args:
        predict_function: The prediction function to call (CDG.predict since it's the first called).
        *args: Arguments to pass to the prediction function.
        max_attempts (int): Maximum number of attempts to make.
        sleep_time (int or float): Time to wait between attempts, in seconds.
    Returns:
        The result of the prediction function if successful.
    Raises:
        Exception: Re-raises the last exception if all attempts fail.
    """
    for attempt in range(max_attempts):
        try:
            # Try to make the prediction
            return predict_function(*args)
        except RuntimeError as e:
            print(f"Attempt {attempt + 1} failed with error: {e}")
            if attempt + 1 < max_attempts:
                print(f"Retrying in {sleep_time} seconds...")
                time.sleep(sleep_time)  # Wait a bit before retrying
            else:
                print("All attempts failed. Raising the last exception.")
                raise  # Re-raise the last exception if out of attempts


def main(models_dir, multifasta_file, output_tsv, output_fasta):
    models_dir = Path(models_dir)

    class_model_dir = models_dir / "classification"
    class_model_path = class_model_dir / "model.p"
    class_vocab_path = class_model_dir / "vocab.pkl"

    annot_model_dir = models_dir / "annotation"
    annot_model_path = annot_model_dir / "model.p"
    annot_vocab_path = annot_model_dir / "vocab.pkl"

    sequences = []

    # Read sequences from the multifasta file
    for record in SeqIO.parse(multifasta_file, "fasta"):
        sequences.append({"sequence": str(record.seq), "name": record.id})

    # Predict class and cleavage for each sequence
    try:
        class_predictions = robust_predict(
            CDG.predict, class_model_path, class_vocab_path, sequences
        )
    except Exception as final_error:
        print(f"Failed to predict class after several attempts: {final_error}")
    cleavage_predictions = ADG.predict(annot_model_path, annot_vocab_path, sequences)

    # The output of nlpprecursor predictions are in JSON format.
    # The code below parses the JSON into TSV and FASTA format.

    fasta_records = []

    with open(output_tsv, "w", newline="\n") as file:
        writer = csv.writer(file, delimiter="\t")

        writer.writerow(
            [
                "name",
                "class",
                "class_score",
                "cleavage_sequence",
                "cleavage_start",
                "cleavage_stop",
                "cleavage_score",
            ]
        )

        for ind, sequence in enumerate(sequences):
            name = sequence["name"]
            class_pred = class_predictions[ind]["class_predictions"][0]
            cleavage_pred = cleavage_predictions[ind]["cleavage_prediction"]

            writer.writerow(
                [
                    name,
                    class_pred["class"],
                    class_pred["score"],
                    cleavage_pred["sequence"],
                    cleavage_pred["start"],
                    cleavage_pred["stop"],
                    cleavage_pred["score"],
                ]
            )


            peptide_id = f"{name}_{class_pred['class']}_{cleavage_pred['start']}_{cleavage_pred['stop']}_nlpprecursor"
            seq_record = SeqRecord(Seq(cleavage_pred["sequence"]), id=peptide_id, description="")
            fasta_records.append(seq_record)

    SeqIO.write(fasta_records, output_fasta, "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run NLPprecursor prediction and output results.')
    parser.add_argument('models_dir', type=str, help='Directory containing model files.')
    parser.add_argument('multifasta_file', type=str, help='Path to input protein multiFASTA file.')
    parser.add_argument('output_tsv', type=str, help='Path to output TSV file.')
    parser.add_argument('output_fasta', type=str, help='Path to output peptide multiFASTA file.')

    args = parser.parse_args()

    main(args.models_dir, args.multifasta_file, args.output_tsv, args.output_fasta)
