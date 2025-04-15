


source("/g/schwab/Marco/repos/trec_metaB/fncts.R")

input_tables <- load_required_tables()

data <- get_main_dataset(input_tables$df_asvs, input_tables$site_mapping, input_tables$mapping) %>%
  mutate(level_0=str_extract(level_1, "Bacteria|Eukaryota|Archaea")) 

tot_reads_per_site <- data %>%
  group_by(site) %>%
  summarize(reads_site=sum(tot_reads))

tot_euk_reads_per_site <- data %>%
  group_by(site, level_0) %>%
  summarize(reads_euk_site=sum(tot_reads))



filtered_one_percent <- data %>%
  left_join(
    tot_euk_reads_per_site
  ) %>%
  mutate(frac=tot_reads/reads_euk_site) %>%
  filter(#str_detect(level_4, "Dino"),
         frac>0.001)





dinos_larger_oone <- filtered_one_percent %>%
  filter(str_detect(level_4, "Dino"),
         frac>0.01)






### test sequence conservation for same level 9 annotation example:

test <- data %>%
  filter(site=="Reykjavik")

## this one has been annoatted 94 times in reykavik
same_spec <- test %>%
  filter(level_9=="Lewinella_sp.")

## extract all unique sequences
all_seq <- input_tables$df_asvs %>%
  filter(asv_id %in% same_spec$asv_id) %>%
  pull(sequence) %>% unique()

## make a consensus sequence and see how many position dont match in between
source("/g/schwab/Marco/repos/sabeRprobes/R/fncts.R")
cons_seq <- make_consensus_sequence(DNAStringSet(x=all_seq))


















