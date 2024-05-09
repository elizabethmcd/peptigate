import argparse
import csv
import sys
import time
from pathlib import Path

import nlpprecursor
import utils
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


def predict_ripp_sequences(models_dir, protein_fasta_file):
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
    - protein_fasta_file (str): The path to the input FASTA file containing sequences to be
      processed.

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

    # Note this read sequence follows recommendations in the NLPPrecursor docs explicitly,
    # as opposed to using the read_fasta() function.
    for record in SeqIO.parse(protein_fasta_file, "fasta"):
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


def extract_ripp_sequences(
    filtered_predictions,
    protein_fasta_file,
    proteins_output_file,
    protein_peptides_output_file,
    predictions_output_file,
    nucleotide_fasta_file=None,
    nucleotides_output_file=None,
    nucleotide_peptides_output_file=None,
):
    """
    Extracts and writes the sequences and their prediction information to specified TSV and FASTA
    files from the filtered predictions.

    Parameters:
    - filtered_predictions (list of tuples): Filtered sequence predictions to be written out.
      Produced by predict_ripp_sequences().
    - protein_fasta_file (str): The path to a FASTA file containing protein sequences.
      This should be the same file used to make the NLPPrecursor predictions.
    - proteins_output_file (str): The path to the output file where protein sequences that gave rise
      to predicted peptides will be saved. Each sequence is written in amino acid FASTA format with
      its ID as the header. Optional.
    - protein_peptides_output_file (str): The path to the output file where peptide sequences will
      be saved. Peptide sequences written to this file are saved in amino acid FASTA format, with
      headers indicating their source ID, start and end positions in the protein sequence, and that
      NLPPrecursor was the source of the annotation.
    - predictions_output_file (str): Path to the output TSV file where predictions will be saved.
    - nucleotide_fasta_file (str): The path to a FASTA file containing nucleotide sequences.
      This should include sequences for the same genes as used to make the NLPPrecursor predictions.
    - nucleotides_output_file (str): The path to the output file where gene sequences that gave rise
      to predicted peptides will be saved. Each sequence is written in FASTA format with its ID as
      the header. Optional but must be provided if nucleotide_fasta_file is provided.
    - nucleotide_peptides_output_file (str): The path to the output file where peptide sequences
      will be saved. Peptide sequences written to this file are saved in nucleotide FASTA format,
      with headers indicating their source ID, start and end positions in the protein sequence, and
      that NLPPrecursor was the source of the annotation. This means this file only contains the
      nucleotide sequence for the peptide itself. Optional but must be provided if
      nucleotide_fasta_file is provided.
    """
    protein_sequences = utils.read_fasta(protein_fasta_file)
    protein_records = []
    protein_peptide_records = []
    predictions = []

    if nucleotide_fasta_file and nucleotides_output_file and nucleotide_peptides_output_file:
        nucleotide_sequences = utils.read_fasta(nucleotide_fasta_file)
        nucleotide_records = []
        nucleotide_peptide_records = []

    for sequence, class_pred, cleavage_pred in filtered_predictions:
        protein_id = sequence["name"]
        peptide_id = f"{protein_id}_start{cleavage_pred['start']}_end{cleavage_pred['stop']}"

        peptide_metadata = {
            "start": cleavage_pred["start"],
            "end": cleavage_pred["stop"],
            "type": "cleavage",
            "class": class_pred["class"],
            "class_score": class_pred["score"],
            "cleavage_score": cleavage_pred["score"],
            "prediction_tool": "nlpprecursor",
        }

        predictions.append(
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

        protein_sequence = protein_sequences.get(protein_id)
        if protein_sequence:
            protein_records.append(SeqRecord(Seq(protein_sequence), id=protein_id, description=""))
            protein_peptide_sequence = cleavage_pred["sequence"]
            description_fields = [f"{key}:{value}" for key, value in peptide_metadata.items()]
            description = " ".join(description_fields)
            protein_peptide_records.append(
                SeqRecord(
                    Seq(protein_peptide_sequence),
                    id=peptide_id,
                    description=description,
                )
            )

        if nucleotide_fasta_file and nucleotides_output_file and nucleotide_peptides_output_file:
            nucleotide_sequence = nucleotide_sequences.get(protein_id)
            # Note that if transdecoder or similar is not used to predict CDS (protein and
            # nucleotide), nucleotide sequences may not have the same ids as protein sequences and
            # this extraction strategy may fail. When that is the case, we don't output anything,
            # and the nucleotide sequence will not be reported for the protein sequence.
            if nucleotide_sequence:
                nucleotide_records.append(
                    SeqRecord(Seq(nucleotide_sequence), id=protein_id, description="")
                )

                nucleotide_peptide_sequence = nucleotide_sequence[
                    cleavage_pred["start"] * 3 : cleavage_pred["stop"] * 3
                ]
                if utils.verify_translation(
                    nucleotide_peptide_sequence, protein_peptide_sequence, to_stop=True
                ):
                    nucleotide_peptide_records.append(
                        SeqRecord(
                            Seq(nucleotide_peptide_sequence),
                            id=peptide_id,
                            description=description,
                        )
                    )

    SeqIO.write(protein_records, proteins_output_file, "fasta")
    SeqIO.write(protein_peptide_records, protein_peptides_output_file, "fasta")
    if nucleotide_fasta_file and nucleotides_output_file and nucleotide_peptides_output_file:
        SeqIO.write(nucleotide_records, nucleotides_output_file, "fasta")
        SeqIO.write(nucleotide_peptide_records, nucleotide_peptides_output_file, "fasta")

    with open(predictions_output_file, "w", newline="") as predictions_out:
        writer = csv.writer(predictions_out, delimiter="\t")
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
        writer.writerows(predictions)


def main(args):
    filtered_predictions = predict_ripp_sequences(args.models_dir, args.protein_fasta_file)

        extract_ripp_sequences(
            filtered_predictions,
            args.protein_fasta_file,
            args.proteins_output_file,
            args.protein_peptides_output_file,
            args.predictions_output_file,
            args.nucleotide_fasta_file,
            args.nucleotides_output_file,
            args.nucleotide_peptides_output_file,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run NLPprecursor prediction and output results.")
    parser.add_argument(
        "--models_dir", required=True, type=str, help="Directory containing model files."
    )
    parser.add_argument(
        "--protein_fasta_file",
        required=True,
        type=str,
        help="The protein FASTA file input to NLPprecursor.",
    )
    parser.add_argument(
        "--proteins_output_file", required=True, type=str, help="The output file path for proteins."
    )
    parser.add_argument(
        "--protein_peptides_output_file",
        required=True,
        type=str,
        help="The output file path for peptides in amino acid format.",
    )
    parser.add_argument(
        "--predictions_output_file",
        required=True,
        type=str,
        help="The output file path for predictions.",
    )
    parser.add_argument(
        "--nucleotide_fasta_file",
        type=str,
        help="The nucleotide FASTA file for the genes (optional).",
    )
    parser.add_argument(
        "--nucleotides_output_file",
        type=str,
        help="The output file path for nucleotide sequences (optional).",
    )
    parser.add_argument(
        "--nucleotide_peptides_output_file",
        type=str,
        help="The output file path for peptides in nucleotide format (optional).",
    )

    args = parser.parse_args()

    main(args)
