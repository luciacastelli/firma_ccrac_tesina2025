# Shiny app para visualizar degs: https://lcastelli.shinyapps.io/degs/
# Load R packages ----
library(shiny)
library(shinythemes)
library(tidyverse)
library(plotly)
library(DT)
load("example.df")

# Set up UI ----
ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel("Visualización de expresión diferencial"),
  sidebarPanel(
    numericInput("pval", "p-valor", 0.05, min = 0, max = 1, step = 0.01), #pval se manda al server
    numericInput("lfc", "logFC", 2, min = 0, max = 100, step = 0.5), # lfc se manda al server
    verbatimTextOutput("n_significant", placeholder = FALSE),
    radioButtons("data_choice", "Seleccioná los datos para mostrar:", # data_choice se manda al server
                 choices = c("Ejemplo (firma)" = "ejemplo", "Subir .csv" = "csv")
                 ),
    conditionalPanel(
      condition = "input.data_choice == 'csv'",
      tags$small('No incluir caracteres especiales (á, í, !...) y respetar las columnas
                 exactamente: "gen", "logFC", "p", y "sentido" (con valores "up" o "down").'),
      tags$br(""),
      fileInput("file", "Subí tu archivo .csv", accept = ".csv"), #file se manda al server
      actionButton("load", "Cargar y procesar datos"), # load se manda al server
      tags$br("")),
    actionButton("reset", "Refrescar tabla", icon = icon("rotate-right", "font-awesome"))
    ),
  mainPanel(
    # agrego css para agrandar el plot a 80% viewport height (vh)
    tags$style(type = "text/css", "#plot {height: 80vh !important;}"),
    tabsetPanel(
      tabPanel("Tabla", DTOutput("table")), #table se manda del server
      tabPanel("Gráfico", plotlyOutput("plot", width = "100%")), #plot se genera en el server, se manda a ui
      tabPanel("NCBI", uiOutput("ncbi_data", width = "100%")), # ncbi_data se genera en el server
      tabPanel("Más Información",
               h4("Cómo subir datos a la web app:"),
               p(HTML('Para subir datos la tabla debe estar en formato .csv 
               y tener las siguientes columnas, respetando los nombres exactamente:
               <strong>"gen"</strong> (con el nombre de cada gen),
               <strong>"logFC"</strong> (con el logaritmo del "fold change"),
               <strong>"p"</strong> (con el p-valor ajustado),
               y <strong>"sentido"</strong> (con valores "up" o "down").
               No debe haber ningún carácter especial (ej. á, í, !...).
               Los valores numéricos tabulados se muestran con 3 cifras significativas.')),
               br(),
               p(HTML('Por cualquier duda o sugerencia comunicarse con <a href = 
                      "mailto:lucicastelli@uade.edu.ar">lucicastelli@uade.edu.ar</a>.'))
               )
      )
    )
  )

# Set up server ----
server <- function(input, output, session) {
  
  # Initialize reactive values for dataset, pval, lfc, and selected gene
  data <- reactiveVal(example.df) # dataset
  pval <- reactive({input$pval})
  lfc <- reactive({input$lfc})
  selected_gene <- reactiveVal(NULL)
  
  # Reactive expression to read the uploaded CSV
  observeEvent(input$load, {
    req(input$file)
    if (input$data_choice == "csv") {
      data(read.csv(input$file$datapath))
      selected_gene(NULL)
    } 
  })
  
  # Volver a cargar toda la tabla de ejemplo al seleccionar la opción
  observeEvent(input$data_choice, {
    if (input$data_choice == "ejemplo") {
      data(example.df)
      selected_gene(NULL)
    } 
  })
  
  # Borrar selección con botón
  observeEvent(input$reset, {
    selected_gene(NULL)
  })
  
  # Update selected gene when plot is clicked
  observeEvent(event_data("plotly_click"), {
    selection <- event_data("plotly_click")
    if (!is.null(selection) && !is.null(selection$customdata)) {
      selected_gene(selection$customdata)
    }
  })
  
  # Output text with number of significant genes
  output$n_significant <- renderText({
    df <- data()
    degs <- df %>% dplyr::filter(p <= pval() & abs(logFC) >= lfc())
    n_up <- sum(degs$sentido == "up")
    n_down <- sum(degs$sentido == "down")
    paste0("Up: ", n_up, ", down: ", n_down)
    })
    
  # Output table
  output$table <- renderDT({
    df <- data()
    if (!is.null(selected_gene())) {
      df <- df[df$gen == selected_gene(),]
    } else {
      # df <- df # reemplazar la línea de abajo si no querés filtración interactiva
      df <- df %>% dplyr::filter(p <= pval() & abs(logFC) >= lfc())
    }
    datatable(df, options = list(pageLength = 10), selection = "single") %>%
      formatSignif(columns = c("logFC", "p"), digits = 3)
    })

  # Output complete plot
  output$plot <- renderPlotly({
    df <- data()
    degs <- df %>% dplyr::filter(p <= pval() & abs(logFC) >= lfc())
    degnames <- degs$gen
    
    volcano <- ggplot(df) +
      aes(y = -log10(p), x = logFC, customdata = gen) +
      geom_point(size = 2, color = ifelse(!df$gen %in% degnames, "grey", ifelse(df$logFC > 0, "red3", "blue3"))) +
      geom_hline(yintercept = -log10(input$pval), linetype = "longdash", colour = "grey", linewidth = 1) +
      geom_vline(xintercept = input$lfc, linetype = "longdash", colour = "#BE684D", linewidth = 1) +
      geom_vline(xintercept = -input$lfc, linetype = "longdash", colour = "#2C467A", linewidth = 1) +
      theme_bw()
    
    selection <- event_data("plotly_click")
    
    if (!is.null(selection)) {
      selected_gene <- selection$customdata
      volcano <- volcano +
        geom_point(size = 2, color = ifelse(df$gen == selected_gene, "green", "grey"))
    } else {
      selected_gene <- NULL
    }
      return(volcano)
  })
  
  # Output NCBI data
  output$ncbi_data <- renderUI({
    if(is.null(selected_gene())) {
      tags$h4("Haz click en un gen del gráfico para ver su información en NCBI. 
              Puede tardar unos segundos.")
    } else {
      ncbi_url <- paste0("https://www.ncbi.nlm.nih.gov/gene/?term=", selected_gene())
      tags$iframe(src = ncbi_url, width = "100%", height = "600px")
    }
  })
   
}

# ShinyApp----
shinyApp(ui = ui, server = server)
