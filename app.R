library(shiny)
library(shinydashboard)
library(tidyverse)
library(plotly)
library(Biostrings)
library(DECIPHER)
library(grid)
library(readxl)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(grid)  # for grobs
library(gridExtra)  # just in case
library(maps)
library(stringdist)


source("/g/schwab/marco/repos/trec_metaB/fncts.R")

# Load data
df_asv_per_sample <- read_tsv("/g/schwab/marco/projects/trec_metaB/prepared_input/asvs_per_sample.tsv")
df_asv_taxonomy <- read_tsv("/g/schwab/marco/projects/trec_metaB/prepared_input/asv_taxonomy.tsv")
df_sampling_sites <- read_tsv("/g/schwab/marco/projects/trec_metaB/prepared_input/sampling_sites.tsv")
df_filtered <- df_asv_taxonomy
# Add count column (1 per sample)
df_sampling_sites <- df_sampling_sites %>%
  mutate(count = 1)

global_env <- list()

ui <- dashboardPage(
  dashboardHeader(title = "TREC metaB Browser"),
  
  dashboardSidebar(
    sidebarMenu(id = "tabs",  # important!
                #menuItem("Overview", tabName = "page1"),
                menuItem("text based search", tabName = "page3"),
                menuItem("sequence based search", tabName = "page4"),
                menuItem("Main analysis", tabName = "page2")
    )
  ),
  
  dashboardBody(
     tabItems(
    #   tabItem(tabName = "page1",
    #           h2("This is Page 1"),
    #           actionButton("go_to_page2", "Go to Page 2 with Plot A")
    #   ),
       tabItem(tabName = "page3",
            h2("Search for Species/Family/Group etc. in the TREC metaB data"),
            
            textInput("text_input", "Enter your input text:", placeholder = "Type something here..."),

            actionButton("text_submit_button", "Search"),
            actionButton("run_analysis_button", "Analyse metaB data"),
            actionButton("go_page3", "Analyze selected taxon"),
            verbatimTextOutput("text_result_output"),  # for showing output of the compilation
            
            ui <- fluidPage(
              titlePanel("Sampling Sites Interactive Plot (Click a Sample)"),
              
                mainPanel(
                  plotlyOutput("treePlot"),
                  plotOutput("treePlot_sec"),
                  plotOutput("main_analysis")
                )
              
            )
            
      ),
      
      tabItem(tabName = "page4",
              h2("Compile something from text input"),
              
              textInput("text_input", "Enter your input text:", placeholder = "Type something here..."),
              
              actionButton("seq_submit_button", "Submit"),
              
              verbatimTextOutput("seq_result_output")  # for showing output of the compilation
      ),
      tabItem(tabName = "page2",
              h2("This is Page 2"),
              ui <- fluidPage(
                actionButton("submit_selection", "get overview for selection"),
                #actionButton("export_main_ov", "Export plot"),
                downloadButton("dl_main_plot", "Download PDF"),
                #downloadButton("export_main_ov", "Download PDF"),
                radioButtons(
                  inputId = "page2_radio_choice",              # The variable name for server access
                  label = "normalize by:",     # The text above the radio buttons
                  choices = c("Eukaryota", "all"), # List of options
                  selected = "Eukaryota",           # Default selected option
                  inline = FALSE                   # TRUE puts them side-by-side
                ),
                radioButtons(
                  inputId = "page2_ylab",              # The variable name for server access
                  label = "Y axis label:",     # The text above the radio buttons
                  choices = c("asv_id", "level_9", "level_8"), # List of options
                  selected = "asv_id",           # Default selected option
                  inline = FALSE                   # TRUE puts them side-by-side
                ),
                radioButtons(
                  inputId = "page2_group_by",              # The variable name for server access
                  label = "Group y-axis by:",     # The text above the radio buttons
                  choices = c("asv_id", "level_9", "level_8"), # List of options
                  selected = "asv_id",           # Default selected option
                  inline = FALSE                   # TRUE puts them side-by-side
                ),
                titlePanel("Sampling Sites Interactive Plot (Click a Sample)"),
                sidebarLayout(
                  sidebarPanel(
                    selectInput("color_by", 
                                label = "Color bars by:",
                                choices = c("size_fraction", "time", "TARA", "depth_m"),
                                selected = "size_fraction"),
                    verbatimTextOutput("selected_sample"),
                    selectInput("level1_filter",
                                label = "Filter by Domain (level_1):",
                                choices = unique(df_asv_taxonomy$level_1),
                                selected = "Eukaryota"),
                    
                    selectInput("level2_filter",
                                label = "Filter by Group (level_2):",
                                choices = NULL,
                                selected = NULL,
                                multiple = FALSE),
                    
                    selectInput("level3_filter",
                                label = "Filter by Subgroup (level_3):",
                                choices = NULL,
                                selected = NULL,
                                multiple = FALSE),
                    
                    selectInput("level4_filter",
                                label = "Filter by level_4:",
                                choices = NULL,
                                selected = NULL,
                                multiple = FALSE),
                    
                    selectInput("level5_filter",
                                label = "Filter by level_5:",
                                choices = NULL,
                                selected = NULL,
                                multiple = FALSE),
                    
                    selectInput("level6_filter",
                                label = "Filter by level_6:",
                                choices = NULL,
                                selected = NULL,
                                multiple = FALSE),
                    
                    selectInput("level7_filter",
                                label = "Filter by level_7:",
                                choices = NULL,
                                selected = NULL,
                                multiple = FALSE),
                    
                    selectInput("level8_filter",
                                label = "Filter by level_8:",
                                choices = NULL,
                                selected = NULL,
                                multiple = FALSE)
                    
                    
                    
                  ),
                  mainPanel(
                    plotOutput("mainPlot")#,
                    # br(),
                    # plotlyOutput("barPlot"),
                    # br(),
                    # plotOutput("secondPlot")
                  )
                )
              )
              
              
      
      )
    )
  )
)

server <- function(input, output, session) { 
  pending_level2 <- reactiveVal(NULL)
  selected_node <- reactiveVal(NULL)
  holdon <- reactiveVal("stop")
  # output$barPlot <- renderPlotly({
  #   p <- ggplot(df_sampling_sites) +
  #     geom_col(aes(
  #       x = site,
  #       y = count,
  #       fill = .data[[input$color_by]],
  #       text = paste("Sample ID:", sample_id)
  #     ),
  #     position = position_stack(),
  #     width = 0.8) +
  #     theme_minimal() +
  #     labs(
  #       title = paste("Samples per Site (colored by", input$color_by, ")"),
  #       x = "Site", y = "Sample Count", fill = input$color_by
  #     )
  #   
  #   ggplotly(p, tooltip = "text") %>%
  #     config(displayModeBar = FALSE)  # Optional: Hide modebar
  # })
  
 

  
  
  # observe({
  #   req(input$level1_filter)
  #   
  #   filtered <- df_asv_taxonomy %>%
  #     filter(level_1 == input$level1_filter)
  #   
  #   level2_choices <- unique(filtered$level_2)
  #   
  #   updateSelectInput(session, "level2_filter",
  #                     choices = c("--all--", level2_choices),
  #                     selected = "--all--")  # Forces no selection
  # })
  # 
  # observe({
  #   req(input$level2_filter)
  #   
  #   filtered <- df_asv_taxonomy %>%
  #     filter(level_1 == input$level1_filter,
  #            level_2 == input$level2_filter)
  #   
  #   level3_choices <- unique(filtered$level_3)
  #   
  #   updateSelectInput(session, "level3_filter",
  #                     choices = c("--all--", level3_choices),
  #                     selected = "--all--")  # Forces no selection
  # })
  
  # FROM HERE code to automatically selects correct selection options for lower levels
  observe({
    req(input$level1_filter, holdon())

    filtered <- df_asv_taxonomy %>%
      filter(level_1 == input$level1_filter)

    level2_choices <- unique(filtered$level_2)

    updateSelectInput(session, "level2_filter",
                      choices = c("--all--", level2_choices),
                      selected = "--all--")
  })

  # level_3
  observe({
    req(input$level2_filter, holdon())

    filtered <- df_asv_taxonomy %>%
      filter(level_1 == input$level1_filter) %>%
      { if (input$level2_filter != "--all--") filter(., level_2 == input$level2_filter) else . }

    level3_choices <- unique(filtered$level_3)

    updateSelectInput(session, "level3_filter",
                      choices = c("--all--", level3_choices),
                      selected = "--all--")
  })

  # level_4
  observe({
    req(input$level3_filter, holdon())

    filtered <- df_asv_taxonomy %>%
      filter(level_1 == input$level1_filter) %>%
      { if (input$level2_filter != "--all--") filter(., level_2 == input$level2_filter) else . } %>%
      { if (input$level3_filter != "--all--") filter(., level_3 == input$level3_filter) else . }

    level4_choices <- unique(filtered$level_4)

    updateSelectInput(session, "level4_filter",
                      choices = c("--all--", level4_choices),
                      selected = "--all--")
  })

  # level_5
  observe({
    req(input$level4_filter, holdon())

    filtered <- df_asv_taxonomy %>%
      filter(level_1 == input$level1_filter) %>%
      { if (input$level2_filter != "--all--") filter(., level_2 == input$level2_filter) else . } %>%
      { if (input$level3_filter != "--all--") filter(., level_3 == input$level3_filter) else . } %>%
      { if (input$level4_filter != "--all--") filter(., level_4 == input$level4_filter) else . }

    level5_choices <- unique(filtered$level_5)

    updateSelectInput(session, "level5_filter",
                      choices = c("--all--", level5_choices),
                      selected = "--all--")
  })

  # level_6
  observe({
    req(input$level5_filter, holdon())

    filtered <- df_asv_taxonomy %>%
      filter(level_1 == input$level1_filter) %>%
      { if (input$level2_filter != "--all--") filter(., level_2 == input$level2_filter) else . } %>%
      { if (input$level3_filter != "--all--") filter(., level_3 == input$level3_filter) else . } %>%
      { if (input$level4_filter != "--all--") filter(., level_4 == input$level4_filter) else . } %>%
      { if (input$level5_filter != "--all--") filter(., level_5 == input$level5_filter) else . }

    level6_choices <- unique(filtered$level_6)

    updateSelectInput(session, "level6_filter",
                      choices = c("--all--", level6_choices),
                      selected = "--all--")
  })

  # level_7
  observe({
    req(input$level6_filter, holdon())

    filtered <- df_asv_taxonomy %>%
      filter(level_1 == input$level1_filter) %>%
      { if (input$level2_filter != "--all--") filter(., level_2 == input$level2_filter) else . } %>%
      { if (input$level3_filter != "--all--") filter(., level_3 == input$level3_filter) else . } %>%
      { if (input$level4_filter != "--all--") filter(., level_4 == input$level4_filter) else . } %>%
      { if (input$level5_filter != "--all--") filter(., level_5 == input$level5_filter) else . } %>%
      { if (input$level6_filter != "--all--") filter(., level_6 == input$level6_filter) else . }

    level7_choices <- unique(filtered$level_7)

    updateSelectInput(session, "level7_filter",
                      choices = c("--all--", level7_choices),
                      selected = "--all--")
  })

  # level_8
  observe({
    req(input$level7_filter, holdon())

    filtered <- df_asv_taxonomy %>%
      filter(level_1 == input$level1_filter) %>%
      { if (input$level2_filter != "--all--") filter(., level_2 == input$level2_filter) else . } %>%
      { if (input$level3_filter != "--all--") filter(., level_3 == input$level3_filter) else . } %>%
      { if (input$level4_filter != "--all--") filter(., level_4 == input$level4_filter) else . } %>%
      { if (input$level5_filter != "--all--") filter(., level_5 == input$level5_filter) else . } %>%
      { if (input$level6_filter != "--all--") filter(., level_6 == input$level6_filter) else . } %>%
      { if (input$level7_filter != "--all--") filter(., level_7 == input$level7_filter) else . }

    level8_choices <- unique(filtered$level_8)

    updateSelectInput(session, "level8_filter",
                      choices = c("--all--", level8_choices),
                      selected = "--all--")
  })
  ## page3 text based query
  # tree_plot_data <- reactiveVal(NULL)
  # 
  # # 2. Update it inside observeEvent
  # observeEvent(input$text_submit_button, {
  #   req(input$text_input)
  #   tp <- tree_taxo(input$text_input, df_asv_taxonomy)
  #   tree_plot_data(tp)
  # })
  # 
  # # 3. Render the plot outside of observeEvent
  # output$treePlot <- renderPlot({
  #   req(tree_plot_data())  # ensure there's something to plot
  #   tree_plot_data()
  # })
  
  # Step 1: Create reactiveVal without assigning yet
  
  
  # Step 2: Update it after user submits input
  observeEvent(input$text_submit_button, {
    req(input$text_input)
    selected_node(input$text_input)
  })
  
  # Step 3: Generate plot based on current node
  
  
  # Step 4: Render plotly plot
  output$treePlot <- renderPlotly({
    req(selected_node())
    ggplotly(
      tree_taxo_new(selected_node(), df_asv_taxonomy), 
      tooltip = "label")
  })
  
  
  # output$treePlot_sec <- renderPlotly({
  #   req(selected_node())
  #   ggplotly(
  #     tree_taxo(selected_node(), df_asv_taxonomy), 
  #     tooltip = "label")
  # })
  
  # Step 5: Handle plotly click to update selected node
  observeEvent(event_data("plotly_click"), {
    click <- event_data("plotly_click")
    print("click recorded")
    print(click)
    if (!is.null(click)) {
      clicked_label <- click$key %||% click$text %||% click$customdata
      #print(clicked_label)
      if (!is.null(clicked_label)) {
        selected_node(clicked_label)
        print(clicked_label)
        updateTextInput(session, "text_input", value = clicked_label)
        p_tree <- tree_taxo(clicked_label, df_asv_taxonomy)
        output$treePlot_sec <- renderPlot(p_tree)
      }
    }
  })
  
  
  observeEvent(input$go_page3, {
    #pending_level2("clicked")   # <- your desired value
    req(selected_node())
    holdon(NULL)
    raw <- tibble(
      level=paste("level", 1:9, sep="_"),
      fallback="--all--"
    )
    
    inheritance <- prepare_tree_taxo_data(selected_node(), df_asv_taxonomy)
    
    last_level <- inheritance %>%
      filter(taxo==selected_node()) %>%
      pull(num_level)
    
    rel_data <- inheritance %>%
      filter(num_level<=last_level)
    
    final_selection_update <- raw %>%
      left_join(rel_data, by="level") %>%
      mutate(
        selection=case_when(is.na(taxo)~fallback,
                            TRUE ~ taxo),
        selection_panel_name=str_replace(level, "_", "") %>% paste0(.,"_filter")
      ) %>%
      select(selection, selection_panel_name, level)
    
    pending_level2(final_selection_update)
    
    updateTabsetPanel(session, "tabs", selected = "page2")
  })
  
  observeEvent(input$tabs, {
    req(input$tabs == "page2", pending_level2())
    
    final_selection_update <- pending_level2()
    
    for (i in 1:9){
      print(final_selection_update[i,])
      updateSelectInput(session,
                        final_selection_update[i,]$selection_panel_name,
                        choices = c("--all--", unique(df_asv_taxonomy[[final_selection_update[i,]$level]])),
                        selected = final_selection_update[i,]$selection)
    }
    print(pending_level2())
    
    pending_level2(NULL)
    #holdon("done")
  }, ignoreInit = TRUE)
  
  ## Step 6: when analyze button is clicked - do analysis based on current node
  observeEvent(input$run_analysis_button, {
    req(input$text_input)
    query <- input$text_input
    rel_colname <- get_highest_score_colname(query, df_asv_taxonomy)
    
    rel_asvs <- df_asv_taxonomy %>%
      select(asv_id, rel_col=all_of(rel_colname)) %>%
      filter(rel_col==query) %>%
      pull(asv_id)
    
    
    asvs_per_sample <- df_asv_per_sample %>%
      filter(asv_id %in% rel_asvs) %>%
      left_join(df_sampling_sites, by="sample_id")
    
    main_abundance_plot <- get_main_abundance_plot(asvs_per_sample, 
                                                   df_asv_taxonomy, 
                                                   df_asv_per_sample,
                                                   norm_type=input$page2_radio_choice,
                                                   ylab=input$page2_ylab)
    
    print(paste("running analysis on:", input$text_input))
    
    output$main_analysis <- renderPlot({
      
      main_abundance_plot
      
    })
    
  })
  
  
  
  
  ## page 4 sequence based query
  observeEvent(input$seq_submit_button, {
    req(input$text_input)  # Ensure it's not empty
    
    # Simulate "compilation" or processing
    result <- paste("You submitted:", input$text_input)
    
    hits <- sequence_query(query_seq = input$text_input, 
                               df_asv_taxonomy)
    
    output$seq_result_output <- renderText({
      paste(length(hits), "hits detected")
    })
  })
  
  main_plot <- reactiveVal(NULL)
  main_plot_nrows <- reactiveVal(NULL)
  main_plot_nsamples <- reactiveVal(NULL)
 
  observeEvent(input$submit_selection, {

    print("SUBMISSIONS")
    
    level_inputs <- c(
      input$level1_filter,
      input$level2_filter,
      input$level3_filter,
      input$level4_filter,
      input$level5_filter,
      input$level6_filter,
      input$level7_filter,
      input$level8_filter
    )
    
    #level_inputs <- c("Eukaryota", "--all--", "--all--", "--all--", "--all--", "--all--", "--all--", "Heterocapsa")
    
    df_filtered <- select_taxo_level(level_inputs, df_asv_per_sample)
    
    rel_asvs <- unique(df_filtered$asv_id)
    
    print(rel_asvs)
    
    asvs_per_sample <- df_asv_per_sample %>%
      filter(asv_id %in% rel_asvs) %>%
      left_join(df_sampling_sites, by="sample_id")
    
    
    print("building main plot")
    print(asvs_per_sample)
    p3_raw <- get_main_abundance_plot(asvs_per_sample, 
                                  df_asv_taxonomy, 
                                  df_asv_per_sample,
                                  input$page2_radio_choice,
                                  ylab=input$page2_ylab,
                                  grouping=input$page2_group_by)
    
    output$mainPlot <- renderPlot({
      
      p3_raw[[1]] 
      
    })
    stopifnot(inherits(p3_raw[[1]], "ggplot"))  # helpful during dev
    main_plot(p3_raw[[1]])
    main_plot_nrows(p3_raw[[2]])
    main_plot_nsamples(nrow(df_sampling_sites))
    #holdon("done")
    #global_env$main_plot <- p3
    
  })
  
  output$mainPlot <- renderPlot({
    req(main_plot())
    main_plot()
  })
  
  # observeEvent(input$export_main_ov, {
  #   
  #   #req(global_env$main_plot)
  #   
  #   print("exporting plot to")
  #   
  #   pdf(file = "/home/rheinnec/test_export.pdf", width=10, height=20)
  #   grid.arrange(global_env$main_plot)
  #   dev.off()
  #   
  #   print("export done")
  #   
  # })
  
  
  # observeEvent(input$export_main_ov, {
  #   req(main_plot())
  #   path <- file.path(getwd(), "test_export.pdf")
  #   tryCatch({
  #     ggplot2::ggsave(path, plot = main_plot(), device = grDevices::cairo_pdf,
  #                     width = 15, height = 20, units = "in", limitsize = FALSE)
  #     showNotification(paste("Exported to", path))
  #   }, error = function(e) {
  #     showNotification(paste("Export failed:", e$message), type = "error")
  #   })
  # })
  output$dl_main_plot <- downloadHandler(
    filename = function() paste0("main_plot_", Sys.Date(), ".pdf"),
    contentType = "application/pdf",
    content = function(file) {
      req(main_plot())
      ggplot2::ggsave(
        filename  = file,
        plot      = main_plot(),         # your cowplot/ggplot object
        device    = grDevices::cairo_pdf,# or "pdf"
        width     = main_plot_nsamples()/6, height = main_plot_nrows()/5, units = "in",
        limitsize = FALSE
      )
    }
  )
  
  
  # Print sample_id on click
  output$selected_sample <- renderPrint({
    click <- event_data("plotly_click")
    if (is.null(click)) return("Click on a sample bar to see its sample_id.")
    
    # Try to extract sample_id from tooltip text
    tooltip <- click$customdata
    if (!is.null(tooltip)) return(tooltip)
    
    # OR: recover sample_id using pointNumber and mapping
    idx <- click$pointNumber + 1
    sample_clicked <- df_sampling_sites$sample_id[idx]
    paste("Clicked Sample ID:", sample_clicked)
  })
  
  # Second plot: example (taxonomy for clicked sample)
  # output$secondPlot <- renderPlot({
  #   click <- event_data("plotly_click")
  #   if (is.null(click)) return(NULL)
  #   
  #   idx <- click$pointNumber + 1
  #   selected_sample <- df_sampling_sites$sample_id[idx]
  #   
  #   print(selected_sample)
  #   
  #   level_inputs <- c(
  #     input$level1_filter,
  #     input$level2_filter,
  #     input$level3_filter,
  #     input$level4_filter,
  #     input$level5_filter,
  #     input$level6_filter,
  #     input$level7_filter,
  #     input$level8_filter
  #   )
  #   
  #   df_filtered <- select_taxo_level(level_inputs, df_asv_per_sample)
  #   
  #   
  #   # Example plot: count by size_fraction for selected site
  #   p2 <- ggplot(df_filtered) +
  #     geom_col(aes(y = col_x, x=nreads, fill = col_fill), show.legend=F) +
  #     theme_minimal() +
  #     labs(
  #       title = paste("Size Fractions for Site:", selected_sample),
  #       x = "Size Fraction", y = "Count"
  #     )
  #   p2
  #   
  #   
  #   
  # })
}

shinyApp(ui, server)
