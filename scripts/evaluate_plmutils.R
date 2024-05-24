library(tidyverse)
library(optparse)
library(Biostrings)

# set up parsing ----------------------------------------------------------

option_list <- list(
  make_option(c("--coding-fasta-file"), type = "character",
              help = "Path to the coding sequences FASTA file"),
  make_option(c("--coding-prediction-file"), type = "character",
              help = "Path to the coding prediction CSV file"),
  make_option(c("--noncoding-fasta-file"), type = "character",
              help = "Path to the noncoding sequences FASTA file"),
  make_option(c("--noncoding-prediction-file"), type = "character",
              help = "Path to the noncoding prediction CSV file"),
  make_option(c("--output-file"), type = "character",
              help = "Path to save the output results", default = "output.tsv")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# define functions --------------------------------------------------------

count_sequences_in_fasta <- function(fasta_file) {
  fasta <- readDNAStringSet(fasta_file)
  length(fasta)
}

calculate_accuracy <- function(TP, TN, FP, FN) {
  (TP + TN) / (TP + TN + FP + FN)
}

calculate_precision <- function(TP, FP) {
  TP / (TP + FP)
}

calculate_recall <- function(TP, FN) {
  TP / (TP + FN)
}

calculate_f_score <- function(TP, FP, FN) {
  precision <- calculate_precision(TP, FP)
  recall <- calculate_recall(TP, FN)
  2 * (precision * recall) / (precision + recall)
}

calculate_mcc <- function(TP, TN, FP, FN) {
  numerator <- TP * TN - FP * FN
  denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  if (denominator == 0) return(0)
  numerator / denominator
}

# read in data ------------------------------------------------------------

coding <- read_csv(args$coding_prediction_file) %>%
  mutate(actual = "coding") %>%
  mutate(predicted = ifelse(predicted_label == "negative", "noncoding", "coding"))

noncoding <- read_csv(args$noncoding_prediction_file) %>%
  mutate(actual = "noncoding") %>%
  mutate(predicted = ifelse(predicted_label == "negative", "noncoding", "coding"))

results <- bind_rows(coding, noncoding)

# adjust tp/tn/fp/fn by number of original input sequences ----------------

coding_before_plmutils_translate <- count_sequences_in_fasta(args$coding_fasta_file)
noncoding_before_plmutils_translate <- count_sequences_in_fasta(args$noncoding_fasta_file)
noncoding_after_plmutils_translate <- nrow(noncoding)
coding_after_plmutils_translate <- nrow(coding)

# calculate TP/TN/FP/FN values --------------------------------------------

summary_results <- results %>%
  group_by(predicted_label, predicted, actual) %>%
  tally()

# coding as coding
TP <- summary_results %>%
  filter(predicted_label == "positive", predicted == "coding", actual == "coding") %>%
  pull(n)

# noncoding as noncoding
TN <- summary_results %>%
  filter(predicted_label == "negative", predicted == "noncoding", actual == "noncoding") %>%
  pull(n)
# account for the transcripts that are actually noncoding but were dropped by plmutils bc of no predicted ORF
if(!all.equal((TN + FP), noncoding_before_plmutils_translate)){
  TN <- TN + (noncoding_before_plmutils_translate - noncoding_after_plmutils_translate) 
}

# noncoding as coding
FP <- summary_results %>%
  filter(predicted_label == "positive", predicted == "coding", actual == "noncoding") %>%
  pull(n)

# coding as noncoding
FN <- summary_results %>%
  filter(predicted_label == "negative", predicted == "noncoding", actual == "coding") %>%
  pull(n)
# account for the transcripts that are actually coding but were dropped by plmutils bc of no predicted ORF
if(!all.equal((FN + TP), coding_before_plmutils_translate)){
  FN <- FN + (coding_before_plmutils_translate - coding_after_plmutils_translate) 
}

# calculate performance ---------------------------------------------------

metrics <- tibble(
  accuracy = calculate_accuracy(TP, TN, FP, FN),
  precision = calculate_precision(TP, FP),
  recall = calculate_recall(TP, FN),
  f1_score = calculate_f_score(TP, FP, FN),
  mcc = calculate_mcc(TP, TN, FP, FN)
)

write_tsv(metrics, args$output_file)
