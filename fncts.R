
load_required_tables <- function(){
  
  mapping <- read_tsv("/g/schwab/Chandni/MetaB_AML_40um/v2/trec-aml-ssuv4v5_dada2_counts.tsv.gz")
  
  df_asvs <- read_tsv("/g/schwab/Chandni/MetaB_AML_40um/v2/trec-aml-ssuv4v5_dada2_asvs.tsv.gz")
  
  trec_id_mapping <- read_csv("/g/schwab/Chandni/MetaB_AML_40um/TREC_Plankton_MetaB_TRECAML_40um_dataset_files_export.csv") %>%
    mutate(
      sample=str_split(Name, "_") %>% map_chr(.,1),
      trec_name=str_replace_all(Samples, "SAMEA", "TREC")
    ) %>%
    select(err_id=sample, trec_id=Samples, trec_name) %>%
    unique()
  
  
  input_samples <- read_tsv("/g/schwab/Marco/projects/trec_metaB/trec_metaB_site_table - Sheet1.tsv")%>%
    mutate(trec_id=paste0("SAMEA", Barcode_ID))
  
  
  error_samples <- read_xlsx("/g/schwab/Marco/projects/trec_metaB/ErroneousReadset_ToRemove.xlsx")
  
  site_order <- c("Villefranche","Reykjavik","Roscoff", "Tallinn","Kristineberg", "Bilbao","Porto")
  
  site_mapping <- trec_id_mapping %>%
    left_join(input_samples) %>%
    select(err_id, site=Site) %>%
    mutate(
      site=case_when(
        str_detect(site, "Iceland") ~ "Reykjavik",
        str_detect(site, "Villef") ~ "Villefranche",
        TRUE ~ site
      ) %>% factor(., levels=site_order),
    )
  
  max_cols <- max(stringr::str_count(df_asvs$pr2_dada2_taxonomy, ";"), na.rm=T) + 1
  return(
    list(
      mapping=mapping,
      df_asvs=df_asvs %>%
        separate_wider_delim(pr2_dada2_taxonomy, delim = ";", 
                             names = paste0("level_", 1:max_cols), 
                             too_few = "align_start"),
      input_samples=input_samples,
      site_mapping=site_mapping,
      trec_id_mapping=trec_id_mapping,
      error_samples=error_samples
    )
  )
  
}

get_main_dataset <- function(df_asvs, site_mapping, mapping){
  
 
  
  taxo_anno <- df_asvs 
  
  
  
  
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
  
  return(plot_data_wide)
  
}
 


prepare_stacked_barplot_data <- function(df_asvs, site_mapping, mapping){
  
  
  
  plot_data <- get_main_dataset(df_asvs, site_mapping, mapping) %>%
    pivot_longer(cols = c("level_1"),
                 names_to = "level", values_to = "taxo") %>%
    ungroup()
  
  
  return(plot_data)
  
}

# df_asvs <- input_tables$df_asvs
# site_mapping <- input_tables$site_mapping
# mapping <- input_tables$mapping

get_stacked_barplot_per_site <- function(df_asvs, site_mapping, mapping, SITE){
  
  
 plot_data_all_sites <- prepare_stacked_barplot_data(df_asvs, site_mapping, mapping) %>%
   mutate(taxo=factor(taxo))
 
 
 taxo_levels <- levels(plot_data_all_sites$taxo)
 
 print(taxo_levels)
 
 plot_data <- plot_data_all_sites %>%
   filter(site==SITE)
 
 
 p <- ggplot(
   data=plot_data %>% group_by(taxo, site, .drop = F) %>% summarize(tot_reads_new=sum(tot_reads)) %>% 
     mutate(taxo=factor(taxo, levels=taxo_levels)) %>%
     filter(site==SITE),
   aes(
     y=site, x=tot_reads_new, fill=factor(taxo, levels = taxo_levels)
   )
 )
   #facet_wrap(~level, nrow=3)+
 
 
  
  return(p)
}










