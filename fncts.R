
asv_overview_per_sample <- function(asv_frac_per_sample, dino_high_abundance, tree, y_annotation="level_7", tree_xlim=0.2){
  
  
  annotations <- tibble(
    
    level=c(0:8),
    name=paste("level", c(9:1), sep="_"),
    colors=c("purple4","purple", "deeppink", "red", "orange", "orange4", "brown","brown4", "black"),
    labels=c("species", "genus", "family", "order","class", "subdivision", "division", "supergroup", "domain")
    
  )
  
  
  rel_asvs_full <- asv_frac_per_sample %>%
    filter(asv_id %in% dino_high_abundance$asv_id)
  # First, get the list of all relevant samples and their sites
  sample_site_map <- rel_asvs_full %>%
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
      rel_asvs_full %>% 
        filter(asv_id %in% dino_high_abundance$asv_id) %>%
        select(sample, asv_id, frac, level_4),
      by = c("sample", "asv_id")
    ) %>%
    mutate(frac = replace_na(frac, 0)) %>%
    as_tibble()
  
  
  
  mx_frac <- max(plot_data$frac)
  
  legend_scale <- c(0,0.01,0.05,0.1,0.2,round(mx_frac*20)/20)
   
  
  if(mx_frac>0.2){
    color_vec <- c("black","blue4","blue3","darkcyan","green4","green","yellow4","yellow","pink", "deeppink")
    val_vec <- c(0,
                 0.001/mx_frac,
                 0.0099/mx_frac,
                 0.01/mx_frac,
                 0.0499/mx_frac,
                 0.05/mx_frac,
                 0.0999/mx_frac,
                 0.1/mx_frac,
                 0.1999/mx_frac,
                 0.2/mx_frac,
                 1)
  } else {
    color_vec <- c("black","blue4","blue3","darkcyan","green4","green","yellow4","yellow")
    val_vec <- c(0,
                 0.001/mx_frac,
                 0.0099/mx_frac,
                 0.01/mx_frac,
                 0.0499/mx_frac,
                 0.05/mx_frac,
                 0.0999/mx_frac,
                 0.1/mx_frac,
                 1)
  }
  
  
  site_count <- rel_asvs_full %>% group_by(asv_id) %>% 
    summarise(all_sites=paste(unique(site), collapse = ";"),
              n_sites=length(unique(site)))
  
  asv_sorter <- plot_data %>%
    group_by(asv_id) %>%
    summarize(
      #n_sites=length(unique(site[which(frac>0)])),
      mn=mean(frac),
      md=median(frac),
      sm=sum(frac)
    ) %>%
    left_join(site_count) %>%
    arrange(desc(md))
  
  ## tree plot
  
  #library(patchwork)
  
  # Create tree plot (flip if you want to match heatmap)
  p_tree_raw <- 
    ggtree(tree) 
    #scale_x_continuous(limits=c(0,1.5))
  
  
  tree_data <- p_tree_raw$data
  
  # Get the tips in the order they appear (top to bottom by y position)
  asv_phylo_order <- tree_data %>%
    filter(isTip) %>%
    arrange(desc(y)) %>%  # use `desc(y)` if you want top-down order like in heatmaps
    pull(label)
  
  p_tree <- p_tree_raw + 
    scale_y_discrete(limits = asv_phylo_order)+
      coord_cartesian(xlim = c(0, tree_xlim))
  
  
 
  hm <- ggplot(plot_data %>%
                 mutate(asv_id = factor(asv_id, levels = asv_phylo_order)))+
    geom_tile(
      aes(x=sample, y=asv_id, fill=as.numeric(frac)),
      #show.legend = F
    )+
    facet_grid(~site, scales = "free_x", space = "free")+
    scale_fill_gradientn(
      name="abundance [%]",
      values=val_vec,
      breaks=legend_scale,
      labels=paste(legend_scale*100),
      colors = color_vec
    )+
    scale_y_discrete(name=paste0("ASV [", annotations[which(annotations$name==y_annotation),]$labels, "]"),
                     breaks=dino_high_abundance$asv_id, labels=dino_high_abundance[[y_annotation]])+
    theme_bw()+
    theme(
      legend.position = "bottom",
      axis.text = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    )
  
  
  annotation_data <- dino_high_abundance %>%
    rowwise() %>%
    mutate(n_na_levels = sum(is.na(c_across(starts_with("level_"))))) %>%
    ungroup()
  
  
  site_count_plot <- 
    ggplot(asv_sorter %>% 
             left_join(annotation_data %>% select(asv_id, n_na_levels) %>% unique()) %>%
             mutate(asv_id=factor(asv_id, levels = asv_phylo_order)), 
           aes(y=asv_id, x=n_sites, fill=as.character(n_na_levels))) +
    geom_bar(stat="identity", #fill="lightblue"
             show.legend = T
    ) +
    #facet_grid(level_4~, scales = "free", space = "free")+
    # facet_wrap(~1)+
    #coord_flip() +
    theme_bw() +
    scale_fill_manual(
      name="annotation",
      breaks=annotations$level,
      values=annotations$colors,
      labels=annotations$labels)+
    scale_y_discrete(breaks=dino_high_abundance$asv_id, labels=dino_high_abundance[[y_annotation]],
                     position = "right")+
    theme(
      legend.position="bottom",
      #axis.text.y.right = element_text(),
      #axis.text.y.left = element_blank(),
      # axis.ticks.y = element_blank(),
      # axis.title.y = element_blank()
    )
  

  
  # Now use cowplot::plot_grid to combine the heatmap and the site count plot
  combined_plot <- 
    plot_grid(
      rel_widths = c(0.1, 1,0.25),
      p_tree,
      hm,           # Your original heatmap plot
      site_count_plot,  # The plot with total number of sites per ASV
      align = "h",  # Align plots horizontally
      ncol = 3,      # Arrange them in 2 columns
      axis = "bt"
    )
  
  return(combined_plot)
  
}


construct_phylogeny <- function(dino_high_abundance, df_asvs, method="nj"){
  
  
  
  all_cons_seqs <- dino_high_abundance %>%
    select(asv_id, level_6) %>%
    unique() %>%
    left_join(df_asvs %>% select(asv_id, sequence))
  
  
  
  dna_set <- DNAStringSet(all_cons_seqs$sequence)
  names(dna_set) <- all_cons_seqs$asv_id
  
  
  
  aligned_seqs <- AlignSeqs(dna_set)
  
  aln_mat <- as.matrix(aligned_seqs)
  
  
  # Identify columns where the proportion of gaps is below a threshold (e.g., keep positions with <50% gaps)
  keep_cols <- apply(aln_mat, 2, function(col) {
    mean(col == "-") < 0.5
  })
  
  # Subset the alignment
  trimmed_mat <- aln_mat[, keep_cols]
  
  # Convert back to DNAStringSet
  trimmed_seqs <- DNAStringSet(apply(trimmed_mat, 1, paste0, collapse=""))
  names(trimmed_seqs) <- names(aligned_seqs)
  
  if(method=="nj"){
    dm <- DistanceMatrix(trimmed_seqs, includeTerminalGaps=FALSE)
    tree <- nj(dm)
  } else {
    phy <- phyDat(as.matrix(trimmed_seqs), type = "DNA")
    
    # Create starting tree
    dm <- dist.ml(phy)
    treeNJ <- NJ(dm)
    
    # Fit a maximum likelihood model
    fit <- pml(treeNJ, data = phy)
    fit_opt <- optim.pml(fit, model = "GTR", optInv = TRUE, optGamma = TRUE, rearrangement = "stochastic")
    
    tree <- fit_opt$tree
    
  }
  return(tree)
}



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










