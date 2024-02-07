import sys
from nlpprecursor.classification.data import DatasetGenerator as CDG
from nlpprecursor.annotation.data import DatasetGenerator as ADG
from pathlib import Path
import nlpprecursor
from Bio import SeqIO

# This allows for backwards compatibility of the pickled models.
sys.modules["protai"] = nlpprecursor

def main(models_dir, multifasta_file, output_tsv):
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
    class_predictions = CDG.predict(class_model_path, class_vocab_path, sequences)
    cleavage_predictions = ADG.predict(annot_model_path, annot_vocab_path, sequences)

    # The output of nlpprecursor predictions are in JSON format.
    # The code below parses the json into a TSV format, which is easier for downstream tools to parse.
    with open(output_tsv, "w") as f:
        f.write("Name\tClass\tClass Score\tCleavage Sequence\tCleavage Start\tCleavage Stop\tCleavage Score\n")

        for i in range(len(sequences)):
            name = sequences[i]["name"]
            class_pred = class_predictions[i]["class_predictions"][0]
            cleavage_pred = cleavage_predictions[i]["cleavage_prediction"]

            f.write(f"{name}\t{class_pred['class']}\t{class_pred['score']}\t"
                    f"{cleavage_pred['sequence']}\t{cleavage_pred['start']}\t"
                    f"{cleavage_pred['stop']}\t{cleavage_pred['score']}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python nlpprecursor.py <models_dir> <multifasta_file> <output_tsv>")
        sys.exit(1)

    models_dir = sys.argv[1]
    multifasta_file = sys.argv[2]
    output_tsv = sys.argv[3]

    main(models_dir, multifasta_file, output_tsv)
