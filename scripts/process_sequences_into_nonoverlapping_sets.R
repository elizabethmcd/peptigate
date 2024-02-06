library(readr)
library(dplyr)
# library(sourmashconsumr)

# functions ---------------------------------------------------------------

split_train_and_test_data <- function(df, fraction = 0.8, seed = 1){
  # sample without replacement to select fraction of the data set as a training data set
  set.seed(1) # set seed so that the sample command gives the same results each time it is run
  df_annotation <- sort(sample(nrow(df), nrow(df)* fraction))
  train <- df[df_annotation, ]
  test <- df[-df_annotation, ]
  return(list(train_df = train, test_df = test))
}

# read in data sets ------------------------------------------------------

fai_col_names <- c("sequence", "length", "offset", "linebases", "linewidth")
clusters <- read_tsv(snakemake@input[['clusters']], col_names = c("rep", "cluster_member"))
coding_validation <- read_tsv(unlist(snakemake@input[['validation_fai']])[1], col_names = fai_col_names)
noncoding_validation <- read_tsv(unlist(snakemake@input[['validation_fai']])[2], col_names = fai_col_names)
coding_traintest <- read_tsv(unlist(snakemake@input[['traintest_fai']])[1], col_names = fai_col_names)
noncoding_traintest <- read_tsv(unlist(snakemake@input[['traintest_fai']])[2], col_names = fai_col_names)

# check overlap between sets ----------------------------------------------

# upset_df <- from_list_to_upset_df(list(noncoding_validation = noncoding_validation$sequence,
#                                        coding_validation = coding_validation$sequence,
#                                        noncoding_traintest = noncoding_traintest$sequence,
#                                        coding_traintest = coding_traintest$sequence,
#                                        cluster_rep = unique(clusters$rep)))

# filter to non-overlapping sets ------------------------------------------

noncoding_validation_filtered <- noncoding_validation %>%
  filter(sequence %in% clusters$rep) %>%                  # keep sequences that had =<80% homology to other sequences
  filter(!sequence %in% noncoding_traintest$sequence) %>% # remove sequences that are in the noncoding train/test set
  filter(!sequence %in% coding_traintest$sequence)        # remove sequences that ensembl now annotates as coding from validation

coding_validation_filtered <- coding_validation %>% 
  filter(sequence %in% clusters$rep) %>%               # keep sequences that had =<80% homology to other sequences
  filter(!sequence %in% coding_traintest$sequence) %>% # remove seqencues that are in the coding train/test set
  filter(!sequence %in% noncoding_traintest$sequence)  # remove sequences that are in the noncoding train/test set

noncoding_traintest_filtered <- noncoding_traintest %>%
  filter(sequence %in% clusters$rep) # keep sequences that had =<80% homology to other sequences

coding_traintest_filtered <- coding_traintest %>%
  filter(sequence %in% clusters$rep) # keep sequences that had =<80% homology to other sequences

# confirm no unexpected overlap between data sets -----------------------
# upset_df_filtered <- from_list_to_upset_df(list(noncoding_validation = noncoding_validation_filtered$sequence,
#                                                 coding_validation = coding_validation_filtered$sequence,
#                                                 noncoding_traintest = noncoding_traintest_filtered$sequence,
#                                                 coding_traintest = coding_traintest_filtered$sequence,
#                                                 cluster_rep = unique(clusters$rep)))

# subsample the coding traintest set to balance sizes  ------------------

# to avoid biases in classification, it's better to have the same number of sequences in the different classes

# select out all sORFs so these won't be filtered out since they are tricky cases
coding_traintest_filtered_sorfs <- coding_traintest_filtered %>%
  filter(length <= 300)

# sample out larger coding sequences
coding_traintest_filtered_other <- coding_traintest_filtered %>%
  filter(length > 300) %>%
  sample_n(size = nrow(noncoding_traintest_filtered) - nrow(coding_traintest_filtered_sorfs))

# combine to the final set of coding sequences that we'll use for training & testing
coding_traintest_keep <- bind_rows(coding_traintest_filtered_other, coding_traintest_filtered_sorfs)

# split training and testing sets -----------------------------------------

coding_traintest_split <- split_train_and_test_data(coding_traintest_keep)
noncoding_traintest_split <- split_train_and_test_data(noncoding_traintest_filtered)

# write out contigs names for each set ------------------------------------
snakemake_output <- unlist(snakemake@output)
print(snakemake_output)
write.table(coding_traintest_split$train_df$sequence, snakemake_output[1], quote = F, col.names = F, row.names = F)
write.table(coding_traintest_split$test_df$sequence, snakemake_output[2], quote = F, col.names = F, row.names = F)
write.table(coding_validation_filtered$sequence, snakemake_output[3], quote = F, col.names = F, row.names = F)
write.table(noncoding_traintest_split$train_df$sequence, snakemake_output[4], quote = F, col.names = F, row.names = F)
write.table(noncoding_traintest_split$test_df$sequence, snakemake_output[5], quote = F, col.names = F, row.names = F)
write.table(noncoding_validation_filtered$sequence, snakemake_output[6], quote = F, col.names = F, row.names = F)
