library(tidyverse)
library(readxl)

source("/g/schwab/Marco/repos/trec_metaB/fncts.R")

## if not accessing via cluster node or jupyterhub, change paths to local schwab server mount location, e.g. Z:/schwab/
## I recommend acces via jupyterhub R server: https://jupyterhub.embl.de


input_tables <- load_required_tables()

chandni_biobank_table <- read_tsv("/g/schwab/marco/projects/trec_metaB/TREC_STARDUST_Biobank_v1.xlsx - MetaB.tsv")

## prepare tables for Flora and AML people

df_asv_per_sample <- input_tables$mapping %>%
  
  left_join(input_tables$trec_id_mapping,
            by=c("sample"="err_id")) %>%
  
  filter(
    !trec_name %in% input_tables$error_samples$`Ref. collaborateur`
  ) %>%
  mutate(sample_id=str_replace_all(trec_name, "TREC", "")) %>%
  select(asv_id, nreads, sample_id)


df_sampling_sites <- input_tables$trec_id_mapping %>%
  filter(
    !trec_name %in% input_tables$error_samples$`Ref. collaborateur`
  )%>%
  left_join(input_tables$site_mapping) %>%
  mutate(sample_id=str_replace_all(trec_name, "TREC", ""))  %>%
  left_join(
    chandni_biobank_table %>% mutate(sample_id=as.character(`Barcode ID`)),
    #by=c("sample_id"="Barcode ID")
  )%>%
  select(sample_id, genoscope_id=err_id, trec_name, samea_name=trec_id, site,
         date=Date, time=Time, size_fraction=Fraction, TARA, depth_m="Depth (m)")


df_asv_taxonomy <- input_tables$df_asvs %>%
  filter(asv_id %in% df_asv_per_sample$asv_id) %>%
  select(-total, -spread)



write_tsv(df_asv_per_sample, file="/g/schwab/marco/projects/trec_metaB/prepared_input/asvs_per_sample.tsv")

write_tsv(df_asv_taxonomy, file="/g/schwab/marco/projects/trec_metaB/prepared_input/asv_taxonomy.tsv")

write_tsv(df_sampling_sites, file="/g/schwab/marco/projects/trec_metaB/prepared_input/sampling_sites.tsv")



## case 1: get all asvs found in two samples (add/remove for more or less)
SAMPLE_OI <- c("ERR14106000", "ERR14106001") 

all_asvs_for_sample_oi <- input_tables$mapping %>%
  filter(sample %in% SAMPLE_OI) %>%
  left_join(input_tables$df_asvs, by="asv_id")


## case 2: get all sample ids an asv was found in

ASV_OI <- c("76c092e3627725674dda633ea9716795")


all_samples_for_asv_oi <- input_tables$df_asvs %>%
  filter(asv_id %in% ASV_OI) %>%
  left_join(input_tables$mapping, by="asv_id")




## case 3: get per site overview

#SOI <- "Roscoff"


# max_cols <- max(stringr::str_count(input_tables$df_asvs$pr2_dada2_taxonomy, ";"), na.rm=T) + 1
# 
# taxo_anno <- input_tables$df_asvs %>%
#   select(asv_id, pr2_dada2_taxonomy) %>%
#   separate_wider_delim(pr2_dada2_taxonomy, delim = ";", 
#                        names = paste0("level_", 1:max_cols), 
#                        too_few = "align_start")
# 
# 
# 
# 
# plot_data_wide <- input_tables$site_mapping %>%
#   #filter(site==SOI) %>%
#   left_join(input_tables$mapping, by=c("err_id"="sample")) %>%
#   group_by(site, asv_id) %>%
#   summarize(
#     tot_reads=sum(nreads)
#   ) %>%
#   arrange(
#     desc(tot_reads)
#   ) %>%
#   left_join(
#   taxo_anno
#     ) 

plot_data <- prepare_stacked_barplot_data(input_tables$df_asvs, input_tables$site_mapping, input_tables$mapping)

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






