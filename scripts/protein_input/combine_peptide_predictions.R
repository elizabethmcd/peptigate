library(tidyverse)
library(optparse)

option_list <- list(
  make_option(c("--nlpprecursor_path"), type="character",
              default="outputs/demo/cleavage/nlpprecursor/nlpprecursor_predictions.tsv", 
              help="Path to NLPprecursor predictions TSV file."),
  make_option(c("--deeppeptide_path"), type="character",
              default="outputs/demo/cleavage/deeppeptide/deeppeptide_predictions.tsv", 
              help="Path to DeepPeptide predictions TSV file."),
  make_option(c("--sorf_path"), type="character",
              default="outputs/demo/sORF/sorf_peptides.tsv", 
              help="Path to sORF predictions TSV file."),
  make_option(c("--faa_tab_path"), type = "character", 
              default="outputs/demo/predictions/peptides_faa.tsv", 
              help="Path to peptide predictions amino acid sequences TSV file."),
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
#' @param sorf_path Path to the sORF predictions TSV file.
#' @param faa_tab_path Path to amino acid peptide predictions as TSV file.
#' 
#'
#' @return A data frame with peptide predictions merged with various annotations.
combine_peptide_predictions <- function(nlpprecursor_path, 
                                        deeppeptide_path, 
                                        sorf_path,
                                        faa_tab_path) {
  
  nlpprecursor <- read_tsv(nlpprecursor_path) %>%
    select(-nlpprecursor_cleavage_sequence)
  
  deeppeptide <- read_tsv(deeppeptide_path)
  
  sorf <- read_csv(sorf_path) %>%
    select(peptide_id = sequence_id) %>%
    mutate(peptide_type = "sORF", 
           prediction_tool = "less_than_100aa")
  
  peptide_predictions <- bind_rows(nlpprecursor, deeppeptide) %>%
    bind_rows(sorf)
  
  peptide_faa <- read_tsv(faa_tab_path, col_names = c("peptide_id", "protein_sequence"))
  
  peptide_predictions <- left_join(peptide_predictions, peptide_sequences, by = "peptide_id")
}

predictions_df<- combine_peptide_predictions(nlpprecursor_path = args$nlpprecursor_path,
                                             deeppeptide_path = args$deeppeptide_path,
                                             sorf_path = args$sorf_path,
                                             faa_tab_path = args$faa_tab_path)
write_tsv(predictions_df, args$output_predictions_path)
