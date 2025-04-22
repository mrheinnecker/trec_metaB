library(tidyverse)
library(readxl)


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
  filter(frac>0.01)%>%
  rowwise() %>%
  mutate(n_na_levels = sum(is.na(c_across(starts_with("level_"))))) %>%
  ungroup()
  
  
  
#  group_by(asv_id) %>%
#  tally() #%>%
  #filter(n>=2)




# First, get the list of all relevant samples and their sites
sample_site_map <- dino_asv_per_sample %>%
  select(sample, site) %>%
  distinct()

# Then, create the full grid
full_grid <- expand.grid(
  sample = unique(sample_site_map$sample),
  asv_id = unique(dino_high_abundance$asv_id),
  stringsAsFactors = FALSE
) %>%
  left_join(sample_site_map, by = "sample")

# Then join the original data and fill NAs with 0
plot_data <- full_grid %>%
  left_join(
    dino_asv_per_sample %>% 
      filter(asv_id %in% dino_high_abundance$asv_id) %>%
      select(sample, asv_id, frac),
    by = c("sample", "asv_id")
  ) %>%
  mutate(frac = replace_na(frac, 0)) %>%
  as_tibble()



mx_frac <- max(plot_data$frac)

legend_scale <- c(0,0.01,0.05,0.1,round(mx_frac*20)/20)

asv_sorter <- plot_data %>%
  group_by(asv_id) %>%
  summarize(
    n_sites=length(unique(site[which(frac>0)])),
    mn=mean(frac),
    md=median(frac),
    sm=sum(frac)
  ) %>%
  arrange(desc(md))


hm <- ggplot(plot_data %>% mutate(asv_id=factor(asv_id, levels=rev(asv_sorter$asv_id))))+
  geom_tile(
    aes(x=sample, y=asv_id, fill=as.numeric(frac)),
    #show.legend = F
  )+
  facet_grid(~site, scales = "free_x", space = "free")+
  scale_fill_gradientn(
    name="abundance [%]",
    values=c(0,
             0.001/mx_frac,
             0.0099/mx_frac,
             0.01/mx_frac,
             0.0499/mx_frac,
             0.05/mx_frac,
             0.0999/mx_frac,
             0.1/mx_frac,
             1),
    breaks=legend_scale,
    labels=paste(legend_scale*100),
    colors = c("black","blue4","blue3","green4","green","yellow","pink", "deeppink")
  )+
  scale_y_discrete(name="ASV genus",
                   breaks=dino_high_abundance$asv_id, labels=dino_high_abundance$level_8)+
  theme_bw()+
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )




site_count_plot <- 
  ggplot(asv_sorter %>% 
           left_join(dino_high_abundance %>% select(asv_id, n_na_levels) %>% unique()) %>%
           mutate(asv_id=factor(asv_id, levels=rev(asv_sorter$asv_id))), 
         aes(y=asv_id, x=n_sites, fill=as.character(n_na_levels))) +
  geom_bar(stat="identity", #fill="lightblue"
           show.legend = T
           ) +
 # facet_wrap(~1)+
  #coord_flip() +
  theme_bw() +
  scale_fill_manual(
    name="annotation",
    breaks=as.character(c(0:7)),
                      values=c("purple4","purple", "deeppink", "red", "orange", "orange4", "brown", "black"),
                      labels=c("species", "genus", "family", "order", "class", "phylum", "kingdom", "domain"))+
  scale_y_discrete(breaks=dino_high_abundance$asv_id, labels=dino_high_abundance$level_8)+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )

# Now use cowplot::plot_grid to combine the heatmap and the site count plot
combined_plot <- 
  plot_grid(
    rel_widths = c(1,0.25),
  hm,           # Your original heatmap plot
  site_count_plot,  # The plot with total number of sites per ASV
  align = "h",  # Align plots horizontally
  ncol = 2,      # Arrange them in 2 columns
  axis = "bt"
)

pdf(file="/g/schwab/Marco/projects/trec_metaB/dinos_per_sample.pdf", width=12, height=7)
print(combined_plot)
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


















