library(tidyverse)


## if not accessing via cluster node or jupyterhub, change paths to local schwab server mount location, e.g. Z:/schwab/
mapping <- read_tsv("/g/schwab/Chandni/MetaB_AML_40um/trec-aml-ssuv4v5_dada2_counts.tsv")

df_asvs <- read_tsv("/g/schwab/Chandni/MetaB_AML_40um/trec-aml-ssuv4v5_dada2_asvs.tsv/trec-aml-ssuv4v5_dada2_asvs.tsv")


## case 1: get all asvs found in two samples (add/remove for more or less)
SAMPLE_OI <- c("ERR14106000", "ERR14106001") 

all_asvs_for_sample_oi <- mapping %>%
  filter(sample %in% SAMPLE_OI) %>%
  left_join(df_asvs, by="asv_id")


## case 2: get all sample ids an asv was found in

ASV_OI <- c("76c092e3627725674dda633ea9716795")


all_samples_for_asv_oi <- df_asvs %>%
  filter(asv_id %in% ASV_OI) %>%
  left_join(mapping, by="asv_id")







