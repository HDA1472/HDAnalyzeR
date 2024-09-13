library(shiny)
library(bslib)
library(dplyr)
library(HDAnalyzeR)  # Assuming your package is called HDAnalyzeR


# This app is only for basic usage and does not include all the features of HDAnalyzeR
ui <- fluidPage(

  titlePanel("Data Preprocessing App"),

  # Main panel with tabs for different sections
  mainPanel(
    tabsetPanel(
      id = "main_tabs",

      # Tab for importing data
      tabPanel("1. Data", value = "tab_import",
               h3("Import Data and Metadata"),
               p("This app operates exactly like the HDAnalyzeR package, but with a shiny interface."),
               p("The data should contain DAid, Assay and NPX columns."),
               fileInput("datafile", "Upload Data File (csv, tsv, txt, rda, rds, xlsx)",
                         accept = c(".csv", ".tsv", ".txt", ".rda", ".rds", ".xlsx")),

               br(),

               p("The metadata should contain DAid and any metadata columns you need."),
               fileInput("metadatafile", "Upload Metadata File (csv, tsv, txt, rda, rds, xlsx)",
                         accept = c(".csv", ".tsv", ".txt", ".rda", ".rds", ".xlsx"))
      ),

      # Tab for preprocessing data and metadata
      tabPanel("2. Preprocessing", value = "tab_preprocess",
               h3("Preprocess Data"),
               selectInput("keep_cols_data", "Columns to Keep in Data", choices = NULL, multiple = TRUE),
               textInput("filter_plates", "Filter Plates (comma-separated)"),
               textInput("filter_assays", "Filter Assays (comma-separated)"),
               checkboxInput("filter_assay_warning", "Filter Assay Warnings", value = FALSE),
               tableOutput("processed_data_preview"),
               tableOutput("data_head"),
               textOutput("data_dimensions"),
               uiOutput("data_column_selector"),
               uiOutput("data_unique_value_dropdown"),

               br(),

               h3("Preprocess Metadata"),
               selectInput("keep_cols_metadata", "Columns to Keep in Metadata", choices = NULL, multiple = TRUE),
               textInput("cohort", "Cohort to Filter (comma-separated)"),
               tableOutput("processed_metadata_preview"),
               textOutput("metadata_dimensions"),
               uiOutput("metadata_column_selector"),
               uiOutput("metadata_unique_value_dropdown")
      ),

      # Tab for PCA and UMAP analysis
      tabPanel("3. PCA and UMAP", value = "tab_analysis",
               h3("PCA Analysis"),
               p("Missing data are imputed with KNN imputation and 5 neighbours."),
               numericInput("pcs", "Number of PCs to Show", value = 5, min = 1),
               textInput("pca_color", "Color by (Column Name)", value = "Disease"),
               checkboxInput("pca_assay", "Plot Assays instead of samples", value = FALSE),
               textInput("pca_palette", "Palette (e.g., cancers12 or custom named vector)", value = NULL),
               actionButton("pca_run", "Run PCA"),
               div(
                 style = "display: flex; overflow-x: auto; gap: 20px;",  # Flex layout with horizontal scroll and gap between plots

                 # Each plot is wrapped in a div with a fixed width and height
                 div(
                   style = "flex: 0 0 33%;",  # Each plot takes up 33% of the viewport width
                   plotOutput("pca_plot1", width = "600px", height = "450px")
                 ),
                 div(
                   style = "flex: 0 0 33%;",  # Each plot takes up 33% of the viewport width
                   plotOutput("pca_plot2", width = "600px", height = "450px")
                 ),
                 div(
                   style = "flex: 0 0 33%;",  # Each plot takes up 33% of the viewport width
                   plotOutput("pca_plot3", width = "600px", height = "450px")
                 )
               ),

               br(),

               h3("UMAP Analysis"),
               p("Missing data are imputed with KNN imputation and 5 neighbours."),
               textInput("umap_color", "Color by (Column Name)", value = "Disease"),
               checkboxInput("umap_assay", "Plot Assays instead of samples", value = FALSE),
               textInput("umap_palette", "Palette (e.g., cancers12 or custom named vector)", value = NULL),
               actionButton("umap_run", "Run UMAP"),
               plotOutput("umap_plot", width = "600px", height = "450px")
      ),

      # Tab for Differential Expression Analysis
      tabPanel("4. Differential Expression", value = "tab_analysis2",
               h3("Differential Expression Analysis"),
               textInput("de_variable", "Variable for DE Analysis", value = "Disease"),
               textInput("de_case", "Case group for DE Analysis (e.g., 'AML')", value = "AML"),
               textInput("de_controls", "Control group for DE Analysis (comma-separated, e.g. CLL, MYEL)", value = "CLL, MYEL"),
               textInput("de_correct", "Covariates to Correct For (e.g., 'Sex, Age')", value = "Sex, Age"),
               textInput("de_correct_type", "Type of Covariates (comma-separated, e.g. factor, numeric)", value = "factor, numeric"),
               textInput("de_only_female", "Groups Only Female (comma-separated, e.g. Group1, Group2)", value = NULL),
               textInput("de_only_male", "Groups Only Male (comma-separated, e.g. 'Group1, Group2')", value = NULL),
               numericInput("de_pval_lim", "P-value Limit", value = 0.05, min = 0, max = 1, step = 0.05),
               numericInput("de_logfc_lim", "Log Fold Change Limit", value = 0, step = 0.1),
               numericInput("de_top_up_prot", "Top Upregulated Proteins", value = 40, step = 1),
               numericInput("de_top_down_prot", "Top Downregulated Proteins", value = 10, step = 1),
               textInput("de_user_defined_proteins", "User-defined Proteins (comma-separated, e.g. FLT3, EPO)", value = NULL),
               textInput("de_palette", "Palette for DE Analysis (e.g., diff_exp or custom named vector)", value = "diff_exp"),
               actionButton("de_run", "Run Differential Expression Analysis"),
               plotOutput("de_volcano_plot", width = "600px", height = "450px")
        )
    ),
    width = 12
  ),

  theme = bs_theme(preset = "darkly",
                   primary = "#883268",
                   secondary = "#F2F2F2")
)

server <- function(input, output, session) {

  convert_palette <- function(palette_input) {
    if (is.null(palette_input) || palette_input == "") {
      return(NULL)
    }

    if (startsWith(palette_input, "c(") && endsWith(palette_input, ")")) {
      palette_input <- substring(palette_input, 3, nchar(palette_input) - 1)  # Remove 'c(' and ')'

      # Convert to named vector
      tryCatch({
        palette <- eval(parse(text = paste0("c(", palette_input, ")")))
        return(palette)
      }, error = function(e) {
        warning("Error parsing palette: ", e$message)
        return(NULL)
      })
      } else {
      # Predefined palette names
      palette <- palette_input
    }
    return(palette)
  }

  # Reactive to load the data file using HDAnalyzeR::import_df
  data <- reactive({
    req(input$datafile)
    HDAnalyzeR::import_df(input$datafile$datapath)
  })

  # Reactive to load the metadata file using HDAnalyzeR::import_df
  metadata <- reactive({
    req(input$metadatafile)
    HDAnalyzeR::import_df(input$metadatafile$datapath)
  })

  # Dynamically update column selections for data
  observe({
    req(data())
    updateSelectInput(session, "keep_cols_data", choices = colnames(data()), selected = colnames(data()))
  })

  # Dynamically update column selections for metadata
  observe({
    req(metadata())
    updateSelectInput(session, "keep_cols_metadata", choices = colnames(metadata()), selected = colnames(metadata()))
  })

  # Process the data using HDAnalyzeR::clean_data
  processed_data <- reactive({
    req(data())

    filter_plates <- if (!is.null(input$filter_plates) && input$filter_plates != "") {
      strsplit(input$filter_plates, ", ")[[1]]
    } else NULL

    filter_assays <- if (!is.null(input$filter_assays) && input$filter_assays != "") {
      strsplit(input$filter_assays, ", ")[[1]]
    } else NULL

    HDAnalyzeR::clean_data(
      df_in = data(),
      keep_cols = input$keep_cols_data,
      filter_plates = filter_plates,
      filter_assays = filter_assays,
      filter_assay_warning = input$filter_assay_warning
    )
  })

  # Process the metadata using HDAnalyzeR::clean_metadata
  processed_metadata <- reactive({
    req(metadata())

    cohort <- if (!is.null(input$cohort) && input$cohort != "") {
      strsplit(input$cohort, ", ")[[1]]
    } else NULL

    HDAnalyzeR::clean_metadata(
      df_in = metadata(),
      keep_cols = input$keep_cols_metadata,
      cohort = cohort
    )
  })

  # Show the first 3 rows of the processed data
  output$processed_data_preview <- renderTable({
    head(processed_data(), 10)
  })

  # Show the number of rows and columns of the processed data
  output$data_dimensions <- renderText({
    df <- processed_data()
    paste("Rows:", nrow(df), "Columns:", ncol(df))
  })

  # Dynamically create a dropdown for column names
  output$data_column_selector <- renderUI({
    selectInput("selected_column", "Choose a Column:", choices = names(processed_data()))
  })

  output$data_unique_value_dropdown <- renderUI({
    req(input$selected_column)  # Make sure a column is selected

    unique_values <- unique(processed_data()[[input$selected_column]])

    selectInput("selected_value", "Unique Values:", choices = unique_values, multiple = FALSE)
  })

  # Show the first 3 rows of the processed metadata
  output$processed_metadata_preview <- renderTable({
    head(processed_metadata(), 10)
  })

  # Show the number of rows and columns of the processed metadata
  output$metadata_dimensions <- renderText({
    df <- processed_metadata()
    paste("Rows:", nrow(df), "Columns:", ncol(df))
  })

  # Dynamically create a dropdown for column names
  output$metadata_column_selector <- renderUI({
    selectInput("selected_column", "Choose a Column:", choices = names(processed_metadata()))
  })

  output$metadata_unique_value_dropdown <- renderUI({
    req(input$selected_column)  # Make sure a column is selected

    unique_values <- unique(processed_metadata()[[input$selected_column]])

    selectInput("selected_value", "Unique Values:", choices = unique_values, multiple = FALSE)
  })

  # Reactive to trigger PCA analysis
  observeEvent(input$pca_run, {
    req(processed_data())

    # Set metadata to NULL if assay is TRUE
    metadata <- if (input$pca_assay) NULL else processed_metadata()
    palette <- convert_palette(input$pca_palette)

    pca_results <- HDAnalyzeR::do_pca(
      olink_data = processed_data(),
      metadata = metadata,
      pcs = input$pcs,
      color = input$pca_color,
      palette = palette,
      wide = FALSE,       # Fixed argument
      assay = input$pca_assay,
      impute = TRUE,     # Fixed argument
      plots = TRUE,      # Fixed argument
      x = "PC1",         # Default
      y = "PC2",         # Default
      npcs = 4,          # Default
      nproteins = 8,     # Default
      loadings = FALSE,  # Default
      save = FALSE       # Default
    )

    output$pca_plots <- renderUI({
      fluidRow(
        column(4, plotOutput("pca_plot1", width = "600px", height = "450px")),
        column(4, plotOutput("pca_plot2", width = "600px", height = "450px")),
        column(4, plotOutput("pca_plot3", width = "600px", height = "450px"))
      )
    })

    output$pca_plot1 <- renderPlot({
      pca_results$pca_plot
    })

    output$pca_plot2 <- renderPlot({
      pca_results$loadings_plot
    })

    output$pca_plot3 <- renderPlot({
      pca_results$variance_plot
    })
  })

  # Reactive to trigger UMAP analysis
  observeEvent(input$umap_run, {
    req(processed_data())

    # Set metadata to NULL if assay is TRUE
    metadata <- if (input$umap_assay) NULL else processed_metadata()
    palette <- convert_palette(input$umap_palette)

    umap_results <- HDAnalyzeR::do_umap(
      olink_data = processed_data(),
      metadata = metadata,
      color = input$umap_color,
      palette = palette,
      wide = FALSE,       # Fixed argument
      assay = input$umap_assay,
      impute = TRUE,     # Fixed argument
      plots = TRUE,      # Fixed argument
      save = FALSE       # Fixed argument
    )

    output$umap_plot <- renderPlot({
      umap_results$umap_plot
    })
  })

  # Differential Expression Analysis
  observeEvent(input$de_run, {
    req(processed_data())

    # Parse inputs
    variable <- input$de_variable
    case <- input$de_case
    controls <- strsplit(input$de_controls, ",\\s*")[[1]]  # Split comma-separated list into vector
    correct <- if (input$de_correct == "") NULL else strsplit(input$de_correct, ",\\s*")[[1]]
    correct_type <- if (input$de_correct_type == "") NULL else strsplit(input$de_correct_type, ",\\s*")[[1]]
    only_female <- strsplit(input$de_only_female, ",\\s*")[[1]]
    only_male <- strsplit(input$de_only_male, ",\\s*")[[1]]
    pval_lim <- input$de_pval_lim
    logfc_lim <- input$de_logfc_lim
    top_up_prot <- input$de_top_up_prot
    top_down_prot <- input$de_top_down_prot
    palette <- convert_palette(input$de_palette)
    user_defined_proteins <- if (input$de_user_defined_proteins == "") NULL else strsplit(input$de_user_defined_proteins, ",\\s*")[[1]]

    de_results <- HDAnalyzeR::do_limma(
      olink_data = processed_data(),
      metadata = processed_metadata(),
      variable = variable,
      case = case,
      control = controls,
      correct = correct,
      correct_type = correct_type,
      wide = FALSE,
      only_female = only_female,
      only_male = only_male,
      volcano = TRUE,
      pval_lim = pval_lim,
      logfc_lim = logfc_lim,
      top_up_prot = top_up_prot,
      top_down_prot = top_down_prot,
      palette = palette,
      report_nproteins = TRUE,
      user_defined_proteins = user_defined_proteins,
      subtitle = NULL,
      save = FALSE
    )

    output$de_volcano_plot <- renderPlot({
      de_results$volcano_plot  # Adjust based on your actual `do_limma` function
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)
