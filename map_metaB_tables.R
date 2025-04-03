library(tidyverse)
library(readxl)

## if not accessing via cluster node or jupyterhub, change paths to local schwab server mount location, e.g. Z:/schwab/
## I recommend acces via jupyterhub R server: https://jupyterhub.embl.de
mapping <- read_tsv("/g/schwab/Chandni/MetaB_AML_40um/trec-aml-ssuv4v5_dada2_counts.tsv")

df_asvs <- read_tsv("/g/schwab/Chandni/MetaB_AML_40um/trec-aml-ssuv4v5_dada2_asvs.tsv/trec-aml-ssuv4v5_dada2_asvs.tsv")

trec_id_mapping <- read_csv("/g/schwab/Chandni/MetaB_AML_40um/TREC_Plankton_MetaB_TRECAML_40um_dataset_files_export.csv") %>%
  mutate(
    sample=str_split(Name, "_") %>% map_chr(.,1)
  ) %>%
  select(err_id=sample, trec_id=Samples) %>%
  unique()


input_samples <- read_tsv("/g/schwab/Marco/projects/trec_metaB/trec_metaB_site_table - Sheet1.tsv")%>%
  mutate(trec_id=paste0("SAMEA", Barcode_ID))


site_mapping <- trec_id_mapping %>%
  left_join(input_samples) %>%
  select(err_id, site=Site)



# 
# input_samples <- read_xlsx("/g/schwab/Chandni/MetaB_AML_40um/Genoscope/Template_for_Sample_MetaB_ChandniKarel.xlsx", skip=1) %>%
#   
#   select(name, trec_id)


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




## case 3: get per site overview

#SOI <- "Roscoff"


max_cols <- max(stringr::str_count(df_asvs$pr2_dada2_taxonomy, ";"), na.rm=T) + 1

taxo_anno <- df_asvs %>%
  select(asv_id, pr2_dada2_taxonomy) %>%
  separate_wider_delim(pr2_dada2_taxonomy, delim = ";", names = paste0("level_", 1:max_cols), too_few = "align_start")




plot_data <- site_mapping %>%
  #filter(site==SOI) %>%
  left_join(mapping, by=c("err_id"="sample")) %>%
  group_by(site, asv_id) %>%
  summarize(
    tot_reads=sum(nreads)
  ) %>%
  arrange(
    desc(tot_reads)
  ) %>%
  left_join(
  taxo_anno
    ) %>%
  pivot_longer(cols = c("level_1", "level_2"),
               names_to = "level", values_to = "taxo")
  


p <- ggplot(
  data=plot_data,
  aes(
  x=site, y=tot_reads, fill=taxo
  )
)+
  facet_wrap(~level, nrow=3)+
  geom_col()+
  theme_bw()+
  theme(
    axis.text.x=element_text(angle=315, hjust=0)
  )

pdf("/g/schwab/Marco/ov_plot.pdf", width=10, height=10)
p
dev.off()














