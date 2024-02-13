library(readr)
library(dplyr)

# functions ---------------------------------------------------------------

split_train_and_test_data <- function(df, fraction = 0.8, seed = 1){
  # sample without replacement to select fraction of the data set as a training data set
  set.seed(seed) # set seed so that the sample command gives the same results each time it is run
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

fai_col_names <- c("sequence", "length", "offset", "linebases", "linewidth")
clusters <- read_tsv("outputs/models/rnasamba/build/1_homology_reduction/clustered_sequences_cluster.tsv", col_names = c("rep", "cluster_member"))
coding_validation <- read_tsv("inputs/validation/rnachallenge/mRNAs.fa.fai", col_names = fai_col_names)
noncoding_validation <- read_tsv("inputs/validation/rnachallenge/ncRNAs.fa.fai", col_names = fai_col_names)
coding_traintest <- read_tsv("outputs/models/rnasamba/build/2_sequence_sets/traintest/all_coding.fa.fai", col_names = fai_col_names)
noncoding_traintest <- read_tsv("outputs/models/rnasamba/build/2_sequence_sets/traintest/all_noncoding.fa.fai", col_names = fai_col_names)
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

# augment the noncoding traintest set to balance sizes  ------------------

# To avoid biases in classification, it's better to have the same number of sequences in the different classes.
# This could also be done with changing the weights on the activation function of the model building in tensorflow,
# but that would require altering the source code of the RNAsamba tool.
# This approach should produce equivalent results.

num_rows_target <- nrow(coding_traintest_filtered)

noncoding_traintest_filtered <- noncoding_traintest_filtered %>%
  slice_sample(n = num_rows_target, replace = TRUE)

# split training and testing sets -----------------------------------------

coding_traintest_split <- split_train_and_test_data(coding_traintest_filtered)
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
