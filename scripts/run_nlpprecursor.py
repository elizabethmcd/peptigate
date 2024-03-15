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
    Secondly, it catches KeyError 'U' for failed predictions.
    Args:
        predict_function: The prediction function to call (CDG.predict since it's the first called).
        *args: Arguments to pass to the prediction function.
        max_attempts (int): Maximum number of attempts to make.
        sleep_time (int or float): Time to wait between attempts, in seconds.
    Returns:
        The result of the prediction function if successful.
        The result 'undetermined' if KeyError 'U' is raised.
        The result 'failed' in any other circumstances.
    """
    for attempt in range(max_attempts):
        try:
            return predict_function(*args)
        except (RuntimeError, KeyError) as e:
            print(f"Attempt {attempt + 1} failed with error: {e}")
            if isinstance(e, KeyError):
                print("KeyError encountered. Labeling as undetermined and continuing.")
                return "undetermined"  # Specific catch for the KeyError to handle the 'U' error
            if attempt + 1 < max_attempts:
                print(f"Retrying in {sleep_time} seconds...")
                time.sleep(sleep_time)
            else:
                print("All attempts failed. Moving to next protein.")
                return "failed"

def predict_ripp_sequences(models_dir, input_fasta):
    """
    Uses NLPPrecursor to predict the class and cleavage sites of sequences from an input FASTA file.
    It filters out sequences classified as "NONRIPP" since these are the negative class.

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

    #sequences = []
    ripp_sequence_predictions = []  # List to hold successful predictions

    # First, process each sequence individually for class predictions
    for record in SeqIO.parse(input_fasta, "fasta"):
        sequence = {"sequence": str(record.seq), "name": record.id}
        prediction = robust_predict(CDG.predict, class_model_path, class_vocab_path, [sequence])
        if prediction != "undetermined" and prediction != "failed":
            if prediction[0]["class_predictions"][0]["class"] != "NONRIPP":
                ripp_sequence_predictions.append((sequence, prediction[0]["class_predictions"][0]))

    # Now, for each successfully predicted sequence, predict the cleavage
    filtered_predictions = []
    for sequence, class_pred in ripp_sequence_predictions:
        cleavage_prediction = ADG.predict(annot_model_path, annot_vocab_path, [sequence])
        cleavage_pred = cleavage_prediction[0]["cleavage_prediction"]
        filtered_predictions.append((sequence, class_pred, cleavage_pred))

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
