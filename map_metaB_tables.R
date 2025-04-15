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
  separate_wider_delim(pr2_dada2_taxonomy, delim = ";", 
                       names = paste0("level_", 1:max_cols), 
                       too_few = "align_start")




plot_data_wide <- site_mapping %>%
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
    ) 

plot_data <- plot_data_wide %>%
  pivot_longer(cols = c("level_2"),
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



### yannicks artistic metaB figure

tree_colors <- tibble(
  taxo=factor(c("TSAR","Haptista","Cryptista","Archaeplastida","Amorphea","Obazoa","Excavata","CRuMs", "other")),
  color=c("#b25545","#b08f89", "#8c8fc6", "#afb96c", "#959ca5","#959ca5", "#917990ff", "#ac9e71ff", "grey")
)

tree_pattern <- paste(tree_colors$taxo, collapse="|")

site_order <- c("Villefranche","Reykjavik","Roscoff", "Tallinn","Kristineberg", "Bilbao","Porto")

pd_mapped_to_tree <- plot_data_wide %>%
  filter(level_1!="Bacteria") %>%
  select(site, tot_reads, level_2) %>%
  pivot_longer(cols = c("level_2"),
               names_to = "level", values_to = "taxo") %>%
  mutate(
    
    ## as the taxonomy annotations from ampliseq do not perfectly match with the one we have in the tree from:
    ## https://www.cell.com/trends/ecology-evolution/fulltext/S0169-5347%2819%2930257-5#f0005
    ## we have to map them manually:
    taxo_new=case_when(
      str_detect(taxo, tree_pattern) ~ str_extract(taxo, tree_pattern),
      TRUE ~ "other"
      #TRUE ~ taxo
    ),
    
    ## to make it nicer
    site_new=case_when(
      str_detect(site, "Iceland") ~ "Reykjavik",
      str_detect(site, "Villef") ~ "Villefranche",
      TRUE ~ site
    ) %>% factor(., levels=site_order),
    
    ## order correctly
    
    
  )
library(showtext)
font_add_google("Gothic A1", "century")

#p2_stack_v <- 
  ggplot(
    data=pd_mapped_to_tree,
    aes(
      x=site_new,
      y=tot_reads, fill=factor(taxo_new, levels=tree_colors$taxo)
    )
  )+
  #facet_wrap(~level, nrow=3)+
  geom_col(position="stack", show.legend = F)+
  scale_fill_manual(breaks = tree_colors$taxo, values = tree_colors$color)+
  scale_  
  #scale_x_discrete()+  
  theme_bw(
        base_family = "century",
    )+
  theme(
    axis.text.x=element_text(angle=325, hjust=0),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),

    #axis.title=element_blank(),
    panel.grid = element_blank()
  )





pdf("/g/schwab/Marco/stack_v_smaller_fonts.pdf", width=9, height=6)
p2_stack_v
dev.off()




png("/g/schwab/Marco/stack_v.png", width=600, height=400)
p2_stack_v
dev.off()





### these not


p2_fill_v <- 
  ggplot(
  data=pd_mapped_to_tree,
  aes(
    x=factor(site_new, levels=rev(sort(unique(pd_mapped_to_tree$site_new)))), y=tot_reads, fill=factor(taxo_new, levels=tree_colors$taxo)
  )
)+
  #facet_wrap(~level, nrow=3)+
  geom_col(position="fill", show.legend = F)+
  scale_fill_manual(breaks = tree_colors$taxo, values = tree_colors$color)+
  #scale_x_discrete()+  
  theme_bw()+
  theme(
    axis.text.x=element_text(angle=325, hjust=0),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank()
  )



p2_fill_h <- 
  ggplot(
    data=pd_mapped_to_tree,
    aes(
      x=site_new, y=tot_reads, fill=factor(taxo_new, levels=rev(tree_colors$taxo))
    )
  )+
  #facet_wrap(~level, nrow=3)+
  geom_col(position="fill", show.legend = F)+
  scale_fill_manual(breaks = tree_colors$taxo, values = tree_colors$color)+
  #scale_x_discrete()+  
  theme_bw()+
  theme(
    #axis.text.y=element_text(angle=325, hjust=0),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank()
  )+
  coord_flip()

p2_stack_h <- 
  ggplot(
    data=pd_mapped_to_tree,
    aes(
      x=site_new, y=tot_reads, fill=factor(taxo_new, levels=rev(tree_colors$taxo))
    )
  )+
  #facet_wrap(~level, nrow=3)+
  geom_col(position="stack", show.legend = F)+
  scale_fill_manual(breaks = tree_colors$taxo, values = tree_colors$color)+
  #scale_x_discrete()+  
  theme_bw()+
  theme(
    #axis.text.x=element_text(angle=325, hjust=0),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank()
  )+
  coord_flip()




# 
# 
# png("/g/schwab/Marco/fill_h.png", width=1600, height=1200, res=300)
# p2_fill_h
# dev.off()


pdf("/g/schwab/Marco/fill_h.pdf", width=3, height=2)
p2_fill_h
dev.off()




pdf("/g/schwab/Marco/fill_v.pdf", width=3, height=2)
p2_fill_v
dev.off()



pdf("/g/schwab/Marco/stack_h.pdf", width=3, height=2)
p2_stack_h
dev.off()






