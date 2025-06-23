library(shiny)
library(tidyverse)
library(plotly)

# Load data
df_asv_per_sample <- read_tsv("/g/schwab/marco/projects/trec_metaB/prepared_input/asvs_per_sample.tsv")
df_asv_taxonomy <- read_tsv("/g/schwab/marco/projects/trec_metaB/prepared_input/asv_taxonomy.tsv")
df_sampling_sites <- read_tsv("/g/schwab/marco/projects/trec_metaB/prepared_input/sampling_sites.tsv")
df_filtered <- df_asv_taxonomy
# Add count column (1 per sample)
df_sampling_sites <- df_sampling_sites %>%
  mutate(count = 1)

# UI
ui <- fluidPage(
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
      plotlyOutput("barPlot"),
      br(),
      plotOutput("secondPlot")
    )
  )
)

# Server
server <- function(input, output, session) {
  output$barPlot <- renderPlotly({
    p <- ggplot(df_sampling_sites) +
      geom_col(aes(
        x = site,
        y = count,
        fill = .data[[input$color_by]],
        text = paste("Sample ID:", sample_id)
      ),
      position = position_stack(),
      width = 0.8) +
      theme_minimal() +
      labs(
        title = paste("Samples per Site (colored by", input$color_by, ")"),
        x = "Site", y = "Sample Count", fill = input$color_by
      )
    
    ggplotly(p, tooltip = "text") %>%
      config(displayModeBar = FALSE)  # Optional: Hide modebar
  })
  
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
  
  # level_2
  observe({
    req(input$level1_filter)
    
    filtered <- df_asv_taxonomy %>%
      filter(level_1 == input$level1_filter)
    
    level2_choices <- unique(filtered$level_2)
    
    updateSelectInput(session, "level2_filter",
                      choices = c("--all--", level2_choices),
                      selected = "--all--")
  })
  
  # level_3
  observe({
    req(input$level2_filter)
    
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
    req(input$level3_filter)
    
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
    req(input$level4_filter)
    
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
    req(input$level5_filter)
    
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
    req(input$level6_filter)
    
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
    req(input$level7_filter)
    
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
  output$secondPlot <- renderPlot({
    click <- event_data("plotly_click")
    if (is.null(click)) return(NULL)
    
    idx <- click$pointNumber + 1
    selected_sample <- df_sampling_sites$sample_id[idx]
    
    print(selected_sample)
    
    # Dummy plot showing taxa counts (replace with real logic if needed)
    df_filtered <- df_asv_per_sample %>%
      filter(sample_id == selected_sample,
             #level_1==input$level1_filter
             ) %>%
      left_join(df_asv_taxonomy) 
    
    if(input$level2_filter == "--all--"){
      df_filtered <- df_filtered %>%
        filter(level_1==input$level1_filter) %>%
        mutate(
          col_x=level_2,
          col_fill=level_3
        )
    } else {
      df_filtered <- df_filtered %>%
        filter(level_1==input$level1_filter,
               level_2==input$level2_filter)%>%
        mutate(
          col_x=level_3,
          col_fill=level_4
        )
    }
    
    print(input$level1_filter)
    
    # Example plot: count by size_fraction for selected site
    ggplot(df_filtered) +
      geom_col(aes(y = col_x, x=nreads, fill = col_fill), show.legend=F) +
      theme_minimal() +
      labs(
        title = paste("Size Fractions for Site:", selected_sample),
        x = "Size Fraction", y = "Count"
      )
  })
}

# Run the app
shinyApp(ui = ui, server = server)
