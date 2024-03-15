import argparse
import csv
import sys
import time
from pathlib import Path

import nlpprecursor
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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
            return predict_function(*args)
        except RuntimeError as e:
            print(f"Attempt {attempt + 1} failed with error: {e}")
            if attempt + 1 < max_attempts:
                print(f"Retrying in {sleep_time} seconds...")
                time.sleep(sleep_time)  # Wait a bit before retrying
            else:
                print("All attempts failed. Raising the last exception.")
                raise  # Re-raise the last exception if out of attempts


def predict_ripp_sequences(models_dir, input_fasta):
    """
    Uses NLPPrecursor to predict the class and cleavage sites of sequences from an input FASTA file.
    It filters out sequences classified as "NONRIPP" since these are the negative class.
    NLPPrecursor predict functions do not return predictions in the same order as supplied or as
    each other. This function accounts for that.
    Further, since NONRIPP is the negative class, if only performs cleavage prediction on peptides
    that are not predicted to be part of the negative class.

    Parameters:
    - models_dir (str): The directory path where the model files are stored. Available for download
      here: https://github.com/magarveylab/nlpprecursor/releases.
    - input_fasta (str): The path to the input FASTA file containing sequences to be processed.

    Returns:
    - list of tuples: Each tuple contains the sequence information, class prediction, and cleavage
      prediction for sequences not classified as "NONRIPP".
    """
    models_dir = Path(models_dir)

    class_model_dir = models_dir / "classification"
    class_model_path = class_model_dir / "model.p"
    class_vocab_path = class_model_dir / "vocab.pkl"

    annot_model_dir = models_dir / "annotation"
    annot_model_path = annot_model_dir / "model.p"
    annot_vocab_path = annot_model_dir / "vocab.pkl"

    sequences = []

    for record in SeqIO.parse(input_fasta, "fasta"):
        sequences.append({"sequence": str(record.seq), "name": record.id})

    try:
        class_predictions = robust_predict(
            CDG.predict, class_model_path, class_vocab_path, sequences
        )
    except Exception as final_error:
        print(f"Failed to predict class after several attempts: {final_error}")
        return []

    class_predictions_dict = {
        pred["name"]: pred["class_predictions"][0] for pred in class_predictions
    }

    # Filter out NONRIPP sequences based on class predictions before cleavage prediction
    # NONRIPP is the negative class, so these are not RIPP peptides.
    nonripp_sequences = [
        seq for seq in sequences if class_predictions_dict[seq["name"]]["class"] != "NONRIPP"
    ]

    if not nonripp_sequences:
        print("All sequences were classified as NONRIPP.")
        return []

    # Predict cleavage only for non-NONRIPP sequences
    cleavage_predictions = ADG.predict(annot_model_path, annot_vocab_path, nonripp_sequences)
    # Convert cleavage_predictions to a dictionary for easier access
    cleavage_predictions_dict = {
        pred["name"]: pred["cleavage_prediction"] for pred in cleavage_predictions
    }

    filtered_predictions = []
    for sequence in nonripp_sequences:
        name = sequence["name"]
        # We already know each sequence here is not NONRIPP, so directly get the predictions
        class_pred = class_predictions_dict[name]
        if name in cleavage_predictions_dict:
            cleavage_pred = cleavage_predictions_dict[name]
            filtered_predictions.append((sequence, class_pred, cleavage_pred))
        else:
            print(f"No cleavage prediction for sequence {name}. Skipping.")

    return filtered_predictions


def extract_ripp_sequences(filtered_predictions, output_tsv, output_fasta):
    """
    Extracts and writes the sequences and their prediction information to specified TSV and FASTA
    files from the filtered predictions.

    Parameters:
    - filtered_predictions (list of tuples): Filtered sequence predictions to be written out.
      Produced by predict_ripp_sequences.py.
    - output_tsv (str): The path to the output TSV file.
    - output_fasta (str): The path to the output FASTA file.
    """
    fasta_records = []

    with open(output_tsv, "w", newline="\n") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(
            [
                "peptide_id",
                "start",
                "end",
                "peptide_type",
                "peptide_class",
                "prediction_tool",
                "nlpprecursor_class_score",
                "nlpprecursor_cleavage_sequence",
                "nlpprecursor_cleavage_score",
            ]
        )

        for sequence, class_pred, cleavage_pred in filtered_predictions:
            protein_id = sequence["name"]
            peptide_id = f"{protein_id}_start{cleavage_pred['start']}_end{cleavage_pred['stop']}"

            writer.writerow(
                [
                    peptide_id,
                    cleavage_pred["start"],
                    cleavage_pred["stop"],
                    "cleavage",
                    class_pred["class"],
                    "nlpprecursor",
                    class_pred["score"],
                    cleavage_pred["sequence"],
                    cleavage_pred["score"],
                ]
            )

            peptide_metadata = {
                "start": cleavage_pred["start"],
                "end": cleavage_pred["stop"],
                "type": "cleavage",
                "class": class_pred["class"],
                "class_score": class_pred["score"],
                "cleavage_score": cleavage_pred["score"],
                "prediction_tool": "nlpprecursor",
            }
            description_fields = [f"{key}:{value}" for key, value in peptide_metadata.items()]
            seq_record = SeqRecord(
                Seq(cleavage_pred["sequence"]),
                id=peptide_id,
                description=" ".join(description_fields),
            )
            fasta_records.append(seq_record)

    SeqIO.write(fasta_records, output_fasta, "fasta")


def main(models_dir, input_fasta, output_tsv, output_fasta):
    filtered_predictions = predict_ripp_sequences(models_dir, input_fasta)
    extract_ripp_sequences(filtered_predictions, output_tsv, output_fasta)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run NLPprecursor prediction and output results.")
    parser.add_argument("models_dir", type=str, help="Directory containing model files.")
    parser.add_argument("input_fasta", type=str, help="Path to input protein multiFASTA file.")
    parser.add_argument("output_tsv", type=str, help="Path to output TSV file.")
    parser.add_argument("output_fasta", type=str, help="Path to output peptide multiFASTA file.")

    args = parser.parse_args()

    main(args.models_dir, args.input_fasta, args.output_tsv, args.output_fasta)
