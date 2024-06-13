library(tidyverse)
library(optparse)

option_list <- list(
  make_option(c("--nlpprecursor_path"), type="character",
              default="outputs/demo/cleavage/nlpprecursor/nlpprecursor_predictions.tsv", 
              help="Path to NLPprecursor predictions TSV file."),
  make_option(c("--deeppeptide_path"), type="character",
              default="outputs/demo/cleavage/deeppeptide/deeppeptide_predictions.tsv", 
              help="Path to DeepPeptide predictions TSV file."),
  make_option(c("--plmutils_path"), type="character",
              default="outputs/demo/sORF/plmutils/predictions.csv", 
              help="Path to plmutils sORF predictions CSV file."),
  make_option(c("--faa_tab_path"), type = "character", 
              default="outputs/demo/predictions/peptides_faa.tsv", 
              help="Path to peptide predictions amino acid sequences TSV file."),
  make_option(c("--ffn_tab_path"), type = "character", 
              default="outputs/demo/predictions/peptides_ffn.tsv", 
              help="Path to peptide predictions nucleotide sequences TSV file."),
  make_option(c("--output_predictions_path"), type="character",
              default="outputs/demo/predictions/peptide_predictions.tsv", 
              help="Path to peptide predictions TSV file.")
)

args <- parse_args(OptionParser(option_list=option_list))

#' Combine Peptide Predictions
#'
#' This function processes peptide predictions by reading data from
#' various text files, performs renaming and selection operations on columns, and
#' binds multiple data frames into a single data frame with all peptide 
#' predictions.
#'
#' @param nlpprecursor_path Path to the NLPprecursor predictions TSV file.
#' @param deeppeptide_path Path to the DeepPeptide predictions TSV file.
#' @param plmutils_path Path to the plmutils predictions CSV file.
#' @param faa_tab_path Path to amino acid peptide predictions as TSV file.
#' @param ffn_tab_path Path to nucleotide peptide predictions as TSV file.
#' 
#'
#' @return A data frame with peptide predictions merged with various annotations.
combine_peptide_predictions <- function(nlpprecursor_path, 
                                        deeppeptide_path, 
                                        plmutils_path,
                                        faa_tab_path,
                                        ffn_tab_path) {
  
  nlpprecursor <- read_tsv(nlpprecursor_path, 
                           col_types = cols(peptide_id = col_character(),
                                            start = col_double(),
                                            end = col_double(),
                                            peptide_type = col_character(),
                                            peptide_class = col_character(),
                                            prediction_tool = col_character(),
                                            nlpprecursor_class_score = col_double(),
                                            nlpprecursor_cleavage_sequence = col_character(),
                                            nlpprecursor_cleavage_score = col_double())) %>%
    select(-nlpprecursor_cleavage_sequence)
  
  deeppeptide <- read_tsv(deeppeptide_path, 
                          col_types = cols(peptide_id = col_character(),
                                           start = col_double(),
                                           end = col_double(),
                                           peptide_type = col_character(),
                                           peptide_class = col_character(),
                                           prediction_tool = col_character()))
  
  plmutils <- read_csv(plmutils_path,
                       col_types = cols(sequence_id = col_character(),
                                        predicted_probability = col_double(),
                                        predicted_label = col_character())) %>%
    filter(predicted_label == "positive") %>%
    select(peptide_id = sequence_id, plmutils_class_probability = predicted_probability) %>%
    mutate(peptide_type = "sORF", 
           prediction_tool = "plmutils")
  
  peptide_predictions <- bind_rows(nlpprecursor, deeppeptide) %>%
    bind_rows(plmutils)
  
  peptide_faa <- read_tsv(faa_tab_path, col_names = c("peptide_id", "protein_sequence"))
  peptide_ffn <- read_tsv(ffn_tab_path, col_names = c("peptide_id", "nucleotide_sequence"))
  peptide_sequences <- left_join(peptide_faa, peptide_ffn) %>%
    select(peptide_id, protein_sequence, nucleotide_sequence)
  
  peptide_predictions <- left_join(peptide_predictions, peptide_sequences, by = "peptide_id")
}

predictions_df<- combine_peptide_predictions(nlpprecursor_path = args$nlpprecursor_path,
                                             deeppeptide_path = args$deeppeptide_path,
                                             plmutils_path = args$plmutils_path,
                                             faa_tab_path = args$faa_tab_path,
                                             ffn_tab_path = args$ffn_tab_path)
write_tsv(predictions_df, args$output_predictions_path)
