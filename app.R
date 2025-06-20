library(shiny)
library(tidyverse)
library(plotly)

# Load data
df_asv_per_sample <- read_tsv("/g/schwab/marco/projects/trec_metaB/prepared_input/asvs_per_sample.tsv")
df_asv_taxonomy <- read_tsv("/g/schwab/marco/projects/trec_metaB/prepared_input/asv_taxonomy.tsv")
df_sampling_sites <- read_tsv("/g/schwab/marco/projects/trec_metaB/prepared_input/sampling_sites.tsv")

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
      verbatimTextOutput("selected_sample")
    ),
    mainPanel(
      plotlyOutput("barPlot"),
      br(),
      plotOutput("secondPlot")
    )
  )
)

# Server
server <- function(input, output) {
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
      filter(sample_id == selected_sample) %>%
      left_join(df_asv_taxonomy)
    
    # Example plot: count by size_fraction for selected site
    ggplot(df_filtered) +
      geom_col(aes(x = level_1, y=nreads, fill = asv_id)) +
      theme_minimal() +
      labs(
        title = paste("Size Fractions for Site:", selected_sample),
        x = "Size Fraction", y = "Count"
      )
  })
}

# Run the app
shinyApp(ui = ui, server = server)
