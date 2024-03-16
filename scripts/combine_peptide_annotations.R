library(tidyverse)
library(optparse)

# Define the command line arguments
option_list <- list(
  make_option(c("--nlpprecursor_path"), type="character",
              default="outputs/cleavage/nlpprecursor/nlpprecursor_predictions.tsv", 
              help="Path to NLPprecursor predictions TSV file."),
  make_option(c("--deeppeptide_path"), type="character",
              default="outputs/cleavage/deeppeptide/predictions.tsv", 
              help="Path to DeepPeptide predictions TSV file."),
  make_option(c("--plmutils_path"), type="character",
              default="outputs/sORF/plmutils/predictions.csv", 
              help="Path to plmutils sORF predictions CSV file."),
  make_option(c("--autopeptideml_dir"), type="character",
              default="outputs/annotation/autopeptideml/", 
              help="Path to directory containing AutoPeptideML TSV files."),
  make_option(c("--deepsig_path"), type="character",
              default='outputs/annotation/deepsig/deepsig.tsv', 
              help="Path to DeepSig annotations TSV file."),
  make_option(c("--peptipedia_path"), type="character",
              default="outputs/annotation/peptipedia/1_blastp/matches.tsv", 
              help="Path to Peptipedia BLAST matches TSV file."),
  make_option(c("--characteristics_path"), type="character",
              default="outputs/annotation/characteristics/peptide_characteristics.tsv", 
              help="Path to peptide characteristics TSV file."),
  make_option(c("--output_predictions_path"), type="character",
              default="outputs/annotation/peptide_predictions.tsv", 
              help="Path to peptide predictions TSV file."),
  make_option(c("--output_annotations_path"), type="character",
              default="outputs/annotation/peptide_annotations.tsv", 
              help="Path to peptide annotations TSV file.")
)

args <- parse_args(OptionParser(option_list=option_list))

#' Combine Peptide Predictions
#'
#' This function processes peptide predictions by reading data from
#' various TSV files, performs renaming and selection operations on columns, and
#' binds multiple data frames into a single data frame with all peptide 
#' predictions.
#'
#' @param nlpprecursor_path Path to the NLPprecursor predictions TSV file.
#' @param deeppeptide_path Path to the DeepPeptide predictions TSV file.
#' @param plmutils_path Path to the plmutils predictions CSV file.
#'
#' @return A data frame with peptide predictions merged with various annotations.
combine_peptide_predictions <- function(nlpprecursor_path, 
                                        deeppeptide_path, 
                                        plmutils_path) {
  
  nlpprecursor <- read_tsv(nlpprecursor_path) %>%
    select(-nlpprecursor_cleavage_sequence)
  
  deeppeptide <- read_tsv(deeppeptide_path)

  plmutils <- read_csv(plmutils_path) %>%
    filter(predicted_label == "positive") %>%
    select(peptide_id = sequence_id) %>%
    mutate(peptide_type = "sORF", 
           prediction_tool = "plmutils")
  
  peptide_predictions <- bind_rows(nlpprecursor, deeppeptide) %>%
    bind_rows(plmutils)
}

#' Combine Peptide Annotations
#'
#' This function processes peptide annotations by reading data from
#' various TSV files, performs renaming and selection operations on columns, and
#' joining multiple data frames into a single data frame with all peptide 
#' annotations.
#' 
#' @param autopeptideml_dir Path to the directory containing AutoPeptideML TSV files.
#' @param deepsig_path Path to the DeepSig annotations TSV file.
#' @param peptipedia_path Path to the Peptipedia BLAST matches TSV file.
#' @param characteristics_path Path to the peptide characteristics TSV file.
combine_peptide_annotations <- function(autopeptideml_dir, deepsig_path,
                                        peptipedia_path, characteristics_path){
  autopeptideml_files <- Sys.glob(paste0(autopeptideml_dir, "/*tsv"))
  autopeptideml <- map_dfr(autopeptideml_files, read_tsv) %>%
    group_by(peptide_id, sequence) %>% 
    summarize(across(everything(), ~ first(na.omit(.))))
  
  deepsig <- read_tsv(deepsig_path) %>%
    select(-tool, -tmp1, -tmp2)
  
  peptipedia <- read_tsv(peptipedia_path) %>%
    rename(peptide_id = qseqid)  %>%
    rename_with(.cols = -peptide_id, 
                function(x) { paste0("peptipedia_blast_", x)}) %>%
    group_by(peptide_id) %>%
    add_count(peptide_id) %>%
    rename(peptipedia_num_hits = n) %>%
    slice_max(peptipedia_blast_bitscore, n = 1) %>%
    slice_min(peptipedia_blast_evalue, n = 1) %>%
    ungroup()
  
  characteristics <- read_tsv(characteristics_path)
  
  peptide_predictions_with_annotations <- peptide_predictions %>%
    left_join(autopeptideml, by = "peptide_id", relationship = "one-to-one") %>%
    left_join(characteristics, by = "peptide_id", relationship = "one-to-one") %>%
    left_join(peptipedia, by = "peptide_id") %>%
    left_join(deepsig, by = "peptide_id")
  
  return(peptide_predictions_with_annotations)
}
predictions_df<- combine_peptide_predictions(nlpprecursor_path = args$nlpprecursor_path,
                                             deeppeptide_path = args$deeppeptide_path,
                                             plmutils_path = args$plmutils_path)
write_tsv(predictions_df, args$output_predictions_path)


annotations_df<- combine_peptide_annotations(autopeptideml_dir = args$autopeptideml_dir,
                                             deepsig_path = args$deepsig_path,
                                             peptipedia_path = args$peptipedia_path,
                                             characteristics_path = args$characteristics_path)
write_tsv(annotations_df, args$output_anotations_path)
