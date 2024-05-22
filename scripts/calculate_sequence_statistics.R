library(tidyverse)
source("scripts/parse_sequence_information.R")

files <- unlist(snakemake@input)
fai_col_names <- c("sequence", "length", "offset", "linebases", "linewidth")
fai <- files %>%
  set_names() %>%
  map_dfr(read_tsv, col_names = fai_col_names, .id = "sequence_set") %>%
  parse_sequence_information_from_seqkit_fai(dataset_type = TRUE)

fai %>%
  group_by(coding_type, set_name) %>%
  tally() %>%
  arrange(set_name) %>%
  write_tsv(snakemake@output[['set_summary']])

fai %>%
  group_by(coding_type, set_name, length_description) %>%
  tally() %>%
  arrange(set_name) %>%
  write_tsv(snakemake@output[['set_length_summary']])

fai %>%
  group_by(coding_type, set_name, length_description, genome) %>%
  tally() %>%
  arrange(set_name) %>%
  pivot_wider(id_cols = c(coding_type, set_name, length_description), names_from = "genome", values_from = "n") %>%
  write_tsv(snakemake@output[['set_length_genome_summary']])
