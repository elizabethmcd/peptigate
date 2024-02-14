library(tidyverse)

parse_sequence_information_from_seqkit_fai <- function(df, dataset_type = FALSE){
  parsed_df <- df %>%
    separate(sequence, into = c("sequence_id", "sequence_type", "chromosome", "gene", 
                                "gene_biotype", "transcript_biotype", "gene_symbol",
                                "description"), 
             sep = " ", remove = FALSE, extra = "merge", fill = "right") %>%
    separate(chromosome, into = c("chromosome_type", "genome", "chromosome_id", 
                                  "start", "end", "strand"), 
             sep = ":", remove = FALSE, extra = "merge", fill = "right") %>%
    mutate(genome = ifelse(genome == "BDGP6", "BDGP6.46", genome)) %>% # fix some Drosophila names
    mutate(genome = ifelse(sequence %in% c("K02E7.12.1", "F29C4.4.1", "C39B10.6b.1", "F53F4.16.1", "F53F4.17.1"), "WBcel235", genome)) %>% # fix some C. elegans genomes
    mutate(length_description = ifelse(length <= 300, "<= 300nt", "> 300nt"))
  
  if(dataset_type == TRUE){
    parsed_df <- parsed_df %>%
      mutate(sequence_set = gsub(".fa.seqkit.fai", "", basename(sequence_set))) %>%
      separate(sequence_set, into = c("coding_type", "set_name"), sep = "_", remove = FALSE)
  }
  
  return(parsed_df)
}
