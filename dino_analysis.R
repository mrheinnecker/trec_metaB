library(tidyverse)
library(readxl)
library(cowplot)
library(Biostrings)
library(DECIPHER)
library(ape)
library(phangorn)
library(ggtree)

source("/g/schwab/Marco/repos/trec_metaB/fncts.R")

input_tables <- load_required_tables()

#data <- get_main_dataset(input_tables$df_asvs, input_tables$site_mapping, input_tables$mapping) %>%
 # mutate(level_0=str_extract(level_1, "Bacteria|Eukaryota|Archaea")) 

dino_asvs <- input_tables$df_asvs %>%
  filter(str_detect(level_4, "Dino")) %>%
  pull(asv_id) %>% unique()

total_reads_per_sample <- input_tables$mapping %>%
  group_by(sample) %>%
  summarise(
    tot_reads=sum(nreads)
  )

to_remove <- input_tables$trec_id_mapping %>% 
  filter(trec_name %in% input_tables$error_samples$`Ref. collaborateur`) %>%
  pull(err_id)

asv_frac_per_sample <- input_tables$mapping %>%
  ## remove the sample wchich have errors for now..
  filter(!sample %in% to_remove) %>%
  left_join(
    total_reads_per_sample
  ) %>%
  #filter(asv_id %in% dino_asvs) %>%
  mutate(frac=nreads/tot_reads) %>%
  left_join(
    input_tables$df_asvs %>% select(asv_id, spread, starts_with("level"))
  ) %>%
  left_join(input_tables$site_mapping, by=c("sample"="err_id"))


dino_asv_per_sample <- asv_frac_per_sample %>%
  filter(asv_id %in% dino_asvs) 



dino_high_abundance <- dino_asv_per_sample %>%
  ## keep only dinos that have at leat once (in one sample) and abundance higher than 1 %
  filter(frac>0.01)

dino_tree <- construct_phylogeny(dino_high_abundance, input_tables$df_asvs)  
  
dino_plot <- asv_overview_per_sample(asv_frac_per_sample, dino_high_abundance, dino_tree, y_annotation="level_6")





euk_high_abundance <- asv_frac_per_sample %>%
  
  filter(str_detect(level_1, "Eukary"),
         frac>0.02)

euk_tree <- construct_phylogeny(euk_high_abundance, input_tables$df_asvs)  


full_plot <- asv_overview_per_sample(asv_frac_per_sample, euk_high_abundance, euk_tree, y_annotation="level_5", tree_xlim = 0.65)




pdf(file="/g/schwab/Marco/projects/trec_metaB/dinos_per_sample.pdf", width=14, height=7)
print(dino_plot)
dev.off()



pdf(file="/g/schwab/Marco/projects/trec_metaB/euk_per_sample.pdf", width=14, height=12)
print(full_plot)
dev.off()





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


















