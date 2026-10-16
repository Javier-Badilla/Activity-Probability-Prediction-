# ============================================================
# Peptide Activity Predictor - VERSION 4 (Desktop-Style App)
# ============================================================

library(shiny)
library(shinyjs)
library(DT)

options(shiny.maxRequestSize = 100 * 1024^2)

source("peptide_analysis_engine_v1.R")

APP_VERSION <- "5.0.0"
APP_NAME    <- "Peptide Activity Predictor"

APP_CREATOR     <- "Javier Badilla"
APP_INSTITUTION <- "Joint doctoral thesis \u2014 Pontificia Universidad Cat\u00f3lica de Valpara\u00edso (PUCV) & Universidad T\u00e9cnica Federico Santa Mar\u00eda (UTFSM)"
APP_SPONSOR     <- "N\u00facleo de Biotecnolog\u00eda Curauma (NBC), Pontificia Universidad Cat\u00f3lica de Valpara\u00edso"
APP_CONTACT     <- "javier.badilla.e@mail.pucv.cl"


# ============================================================
# BRANDING
# Logo files (utfsm_logo.png, pucv_logo.png, nbc_logo.png) must be
# placed in a "www" subfolder next to this app file - Shiny serves
# anything under www/ automatically at that relative path. Each logo
# was flattened onto a white background so dark text/lines remain
# legible regardless of the app's background color.
# ============================================================

# A simple inline SVG favicon (flask icon), embedded as a data URI so
# no separate www/ folder or .ico file is required to get a branded
# browser-tab icon. Swap this for a real .ico in www/favicon.ico for
# a permanent icon once you have proper artwork.
.favicon_svg <- paste0(
  "data:image/svg+xml,",
  utils::URLencode(paste0(
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64">',
    '<rect width="64" height="64" rx="12" fill="%234a7c59"/>',
    '<path d="M26 10h12v16l10 22a4 4 0 0 1-4 6H20a4 4 0 0 1-4-6l10-22z" ',
    'fill="none" stroke="white" stroke-width="3" stroke-linejoin="round"/>',
    '<circle cx="32" cy="44" r="3" fill="white"/>',
    '</svg>'
  ), reserved = TRUE)
)

.custom_css <- "
  .navbar-default { background-color: #2f4f4f; border-color: #2f4f4f; }
  .navbar-default .navbar-brand { color: #ffffff !important; font-weight: 600; }
  .navbar-default .navbar-nav > li > a { color: #e8f0e8 !important; }
  .navbar-default .navbar-nav > li > a:hover,
  .navbar-default .navbar-nav > .active > a { background-color: #4a7c59 !important; color: #fff !important; }
  #splash_screen {
    position: fixed; top: 0; left: 0; width: 100%; height: 100%;
    background-color: #2f4f4f; color: white; z-index: 9999;
    display: flex; flex-direction: column; align-items: center; justify-content: center;
  }
  #splash_screen h1 { font-weight: 300; letter-spacing: 1px; }
  #splash_screen .splash-logos { display: flex; gap: 18px; margin-top: 22px; background: white;
    padding: 12px 20px; border-radius: 8px; }
  #splash_screen .splash-logos img { height: 55px; }
  .settings-section { background-color: #f7f9f7; border-radius: 6px; padding: 15px; margin-bottom: 15px; }
  .app-footer { color: #888; font-size: 0.8em; text-align: center; margin-top: 20px; }
  .navbar-logo-strip { display: flex; align-items: center; gap: 10px; padding: 4px 10px;
    background: white; border-radius: 4px; margin: 8px 15px 0 0; float: right; }
  .navbar-logo-strip img { height: 32px; }
  .about-logo-box { background: white; border-radius: 8px; padding: 20px; text-align: center;
    box-shadow: 0 1px 4px rgba(0,0,0,0.1); margin-bottom: 20px; }
  .about-logo-box img { max-height: 90px; max-width: 100%; margin: 10px 20px; }
  .about-credit-box { background-color: #f7f9f7; border-left: 4px solid #4a7c59; border-radius: 4px;
    padding: 15px 20px; margin-top: 15px; }
"


# ============================================================
# SMALL HELPERS (same as Version 3)
# ============================================================

label_with_tooltip <- function(label, tooltip) {
  tags$span(label, title = tooltip,
             style = "border-bottom: 1px dotted #888; cursor: help;")
}

validate_inputs <- function(training_file, new_file, target_class,
                             one_vs_all, comparison_class,
                             n_runs, prob_threshold, auc_lower, corr_threshold) {
  errors <- character(0)

  if (is.null(training_file)) errors <- c(errors, "Please select a training Excel file.")
  if (is.null(new_file))      errors <- c(errors, "Please select a new-peptide Excel file.")

  if (is.null(target_class) || trimws(target_class) == "") {
    errors <- c(errors, "Target class cannot be empty.")
  }
  if (!isTRUE(one_vs_all) && (is.null(comparison_class) || trimws(comparison_class) == "")) {
    errors <- c(errors, "Comparison class cannot be empty (or check 'one-vs-all').")
  }
  if (!isTRUE(one_vs_all) && !is.null(target_class) && !is.null(comparison_class) &&
      trimws(target_class) == trimws(comparison_class)) {
    errors <- c(errors, "Target class and comparison class cannot be the same.")
  }

  if (is.null(n_runs) || is.na(n_runs) || n_runs < 1 || n_runs != round(n_runs)) {
    errors <- c(errors, "Number of stability runs must be a positive whole number.")
  }
  if (is.null(prob_threshold) || is.na(prob_threshold) || prob_threshold <= 0 || prob_threshold >= 1) {
    errors <- c(errors, "Probability threshold must be between 0 and 1 (exclusive).")
  }
  if (is.null(auc_lower) || is.na(auc_lower) || auc_lower <= 0 || auc_lower >= 1) {
    errors <- c(errors, "Minimum feature AUC must be between 0 and 1 (exclusive).")
  }
  if (is.null(corr_threshold) || is.na(corr_threshold) || corr_threshold <= 0 || corr_threshold >= 1) {
    errors <- c(errors, "Correlation threshold must be between 0 and 1 (exclusive).")
  }

  errors
}


# ============================================================
# UI
# ============================================================

ui <- navbarPage(

  title = paste0(APP_NAME, "  v", APP_VERSION),
  id = "main_nav",
  windowTitle = APP_NAME,

  header = tagList(
    useShinyjs(),
    tags$head(
      tags$link(rel = "icon", href = .favicon_svg),
      tags$style(.custom_css)
    ),

    # ---- Startup splash screen (fades out shortly after the app is ready) ----
    tags$div(
      id = "splash_screen",
      tags$h1(APP_NAME),
      tags$p(paste0("Version ", APP_VERSION, " - Loading...")),
      tags$div(class = "spinner-border text-light", role = "status"),
      tags$div(class = "splash-logos",
        tags$img(src = "utfsm_logo.png"),
        tags$img(src = "pucv_logo.png"),
        tags$img(src = "nbc_logo.png")
      )
    ),

    # ---- Small logo strip shown in the navbar on every tab ----
    tags$div(class = "navbar-logo-strip",
      tags$img(src = "nbc_logo.png"),
      tags$img(src = "pucv_logo.png"),
      tags$img(src = "utfsm_logo.png")
    )
  ),

  # ================= TAB 1: RUN =================
  tabPanel("Run Analysis",
    br(),
    fluidRow(
      column(6,
        h4("1. Select Input Files"),
        fileInput("training_file", "Training Excel file (ID | Type | descriptors)",
                  accept = c(".xlsx", ".xls")),
        fileInput("new_file", "New-peptide Excel file (ID | descriptors)",
                  accept = c(".xlsx", ".xls"))
      ),
      column(6,
        h4("2. Run"),
        p("Parameters are configured on the Settings tab. Current target: ",
          strong(textOutput("current_target_display", inline = TRUE)),
          " vs. ", strong(textOutput("current_comparison_display", inline = TRUE))),
        fluidRow(
          column(6, actionButton("run_btn", "Run Analysis", class = "btn-primary", width = "100%")),
          column(6, actionButton("reset_btn", "Reset", width = "100%"))
        )
      )
    ),
    hr(),
    h4("Status"),
    verbatimTextOutput("status"),
    uiOutput("validation_errors"),
    br(),
    p(em("See the Results tab once a run completes."))
  ),

  # ================= TAB 2: SETTINGS =================
  tabPanel("Settings",
    br(),
    div(class = "settings-section",
      h4("Analysis Parameters"),

      div(label_with_tooltip("Target class",
          "The class label in the 'Type' column you want the model to detect (e.g. AFP).")),
      textInput("target_class", NULL, value = "AFP"),

      checkboxInput("one_vs_all", "Compare target vs. ALL other classes (one-vs-all)",
                    value = FALSE),
      conditionalPanel(
        condition = "!input.one_vs_all",
        div(label_with_tooltip("Comparison class",
            "The class label to compare the target against (e.g. FUNG). Ignored if one-vs-all is checked.")),
        textInput("comparison_class", NULL, value = "FUNG")
      ),

      div(label_with_tooltip("Number of stability runs",
          "How many times the full RF/NN/NN training loop repeats. Higher = more stable estimates, but proportionally longer runtime.")),
      numericInput("n_runs", NULL, value = 20, min = 1, step = 1),

      div(label_with_tooltip("Probability threshold",
          "Probability cutoff used for the 'frequency above threshold' stability statistic.")),
      numericInput("prob_threshold", NULL, value = 0.90, min = 0.01, max = 0.99, step = 0.01),

      div(label_with_tooltip("Minimum feature AUC",
          "Descriptors with a per-feature AUC below this are dropped during feature selection (auto-relaxes to 0.55 if too few features remain).")),
      numericInput("auc_lower", NULL, value = 0.60, min = 0.01, max = 0.99, step = 0.01),

      div(label_with_tooltip("Correlation threshold",
          "Descriptors correlated above this value with another descriptor are filtered out.")),
      numericInput("corr_threshold", NULL, value = 0.90, min = 0.01, max = 0.99, step = 0.01)
    ),
    div(class = "settings-section",
      h4("Output"),
      div(label_with_tooltip("Output folder",
          "Base folder for results. A timestamped subfolder is created automatically for each run so past results are never overwritten.")),
      textInput("output_directory", NULL, value = "output")
    )
  ),

  # ================= TAB 3: RESULTS =================
  tabPanel("Results",
    br(),
    tabsetPanel(
      id = "results_tabs",

      tabPanel("Summary",
        br(),
        h4("Dataset Summary"),
        tableOutput("dataset_summary"),
        h4("Overview"),
        verbatimTextOutput("results_overview")
      ),

      tabPanel("Model Performance",
        br(),
        h4("Mean AUC per model"),
        tableOutput("model_performance"),
        h4("Stability Summary (per new peptide)"),
        DTOutput("stability_table")
      ),

      tabPanel("Predictions",
        br(),
        tabsetPanel(
          tabPanel("Simple",          br(), DTOutput("simple_predictions_table")),
          tabPanel("Detailed",        br(), DTOutput("detailed_predictions_table")),
          tabPanel("High Confidence", br(), DTOutput("high_confidence_table")),
          tabPanel("Consensus",       br(), DTOutput("consensus_table"))
        )
      ),

      tabPanel("Features",
        br(),
        h4("RF Feature Importance (from the last stability run)"),
        DTOutput("importance_table"),
        h4("Top-5 Feature Statistics"),
        tableOutput("feature_stats_table")
      ),

      tabPanel("PCA",
        br(),
        h4("3D PCA - Training Data"),
        uiOutput("pca_train_view"),
        h4("3D PCA - Training + New Samples"),
        uiOutput("pca_combined_view")
      ),

      tabPanel("Files / Output",
        br(),
        h4("Output location"),
        verbatimTextOutput("output_dir_text"),
        fluidRow(
          column(4, uiOutput("download_buttons")),
          column(4, uiOutput("open_folder_button"))
        ),
        hr(),
        fluidRow(
          column(8, h4("Run History (this session)")),
          column(4, actionButton("clear_history_btn", "Clear History", width = "100%"))
        ),
        DTOutput("run_history_table")
      ),

      tabPanel("Log",
        br(),
        tags$div(
          style = "max-height: 500px; overflow-y: auto; background-color: #f7f7f7; padding: 8px;",
          verbatimTextOutput("run_log")
        )
      )
    )
  ),

  # ================= TAB 4: ABOUT =================
  tabPanel("About",
    br(),
    fluidRow(
      column(8, offset = 2,

        div(class = "about-logo-box",
          tags$img(src = "nbc_logo.png"),
          tags$img(src = "pucv_logo.png"),
          tags$img(src = "utfsm_logo.png")
        ),

        h3(APP_NAME),
        p(strong("Version: "), APP_VERSION, " (Desktop-style Shiny application)"),
        p(strong("Engine: "), "peptide_analysis_engine.R"),

        div(class = "about-credit-box",
          p(strong("Created by: "), APP_CREATOR),
          p(strong("Developed as part of: "), APP_INSTITUTION),
          p(strong("Developed for: "), APP_SPONSOR),
          p(strong("Contact: "), APP_CONTACT)
        ),

        hr(),
        h4("What this does"),
        p("Predicts peptide biological activity using an ensemble of three models: ",
          "a Random Forest, a Neural Network trained on PCA-reduced descriptors, and ",
          "a Neural Network trained on ANOVA-selected descriptors. Predictions are ",
          "combined into an AUC-weighted ensemble with Platt-scaled probability ",
          "calibration and multi-run stability analysis."),
        h4("Methodology notes"),
        tags$ul(
          tags$li("Train/test split happens before feature selection to avoid leakage."),
          tags$li("Feature selection uses per-descriptor AUC, with an automatic fallback ",
                  "to a relaxed threshold if too few features pass."),
          tags$li("SMOTE oversampling is applied automatically when class imbalance exceeds 3:1."),
          tags$li("RF feature importance reflects the last of the stability runs, not an average."),
        ),
        tags$div(class = "app-footer",
          paste0(APP_NAME, " v", APP_VERSION, " \u2014 ", APP_CREATOR))
      )
    )
  )
)


# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {

  # ---- Hide splash screen shortly after the app is fully ready ----
  session$onFlushed(function() {
    shinyjs::delay(3000, shinyjs::hide("splash_screen", anim = TRUE, animType = "fade", time = 0.4))
  }, once = TRUE)

  rv <- reactiveValues(
    run_output  = NULL,
    running     = FALSE,
    history     = data.frame(
      Timestamp = character(0), Target = character(0), Comparison = character(0),
      N_Runs = integer(0), RF_AUC = numeric(0), NN_PCA_AUC = numeric(0),
      NN_ANOVA_AUC = numeric(0), High_Confidence = integer(0),
      Output_Folder = character(0), stringsAsFactors = FALSE
    )
  )

  # ---- Live display of current target/comparison on the Run tab ----
  output$current_target_display <- renderText({ trimws(input$target_class) })
  output$current_comparison_display <- renderText({
    if (isTRUE(input$one_vs_all)) "ALL (one-vs-all)" else trimws(input$comparison_class)
  })

  # ---- Reset button ----
  observeEvent(input$reset_btn, {
    updateTextInput(session, "target_class", value = "AFP")
    updateCheckboxInput(session, "one_vs_all", value = FALSE)
    updateTextInput(session, "comparison_class", value = "FUNG")
    updateNumericInput(session, "n_runs", value = 20)
    updateNumericInput(session, "prob_threshold", value = 0.90)
    updateNumericInput(session, "auc_lower", value = 0.60)
    updateNumericInput(session, "corr_threshold", value = 0.90)
    updateTextInput(session, "output_directory", value = "output")
    rv$run_output <- NULL
  })

  # ---- Clear run history ----
  observeEvent(input$clear_history_btn, {
    rv$history <- rv$history[0, ]
  })

  # ---- Run button ----
  observeEvent(input$run_btn, {

    if (isTRUE(rv$running)) return()

    errors <- validate_inputs(
      input$training_file, input$new_file, input$target_class,
      input$one_vs_all, input$comparison_class,
      input$n_runs, input$prob_threshold, input$auc_lower, input$corr_threshold
    )

    if (length(errors) > 0) {
      rv$run_output <- list(success = FALSE, error = NULL, validation_errors = errors, log = "")
      return()
    }

    rv$running <- TRUE
    shinyjs::disable("run_btn")
    shinyjs::disable("reset_btn")
    on.exit({
      rv$running <- FALSE
      shinyjs::enable("run_btn")
      shinyjs::enable("reset_btn")
    })

    run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_dir_this_run <- file.path(trimws(input$output_directory), run_timestamp)

    progress <- shiny::Progress$new(session, min = 0, max = 1)
    progress$set(message = "Running analysis...", value = 0)
    on.exit(progress$close(), add = TRUE)

    progress_callback <- function(value, detail) {
      progress$set(value = value, detail = detail)
    }

    dir.create(output_dir_this_run, recursive = TRUE, showWarnings = FALSE)
    log_file <- file.path(output_dir_this_run, "run_log.txt")
    
    run_output <- tryCatch({
      sink(log_file, split = TRUE)   # writes to file AND to the R console
      analysis_result <- run_analysis(
        training_file     = input$training_file$datapath,
        new_file          = input$new_file$datapath,
        target_class      = trimws(input$target_class),
        comparison_class  = if (isTRUE(input$one_vs_all)) NULL else trimws(input$comparison_class),
        n_runs            = round(input$n_runs),
        prob_threshold    = input$prob_threshold,
        auc_lower         = input$auc_lower,
        corr_threshold    = input$corr_threshold,
        output_directory  = output_dir_this_run,
        progress_callback = progress_callback
      )
      sink()
      log_lines <- readLines(log_file, warn = FALSE)
      list(success = TRUE, result = analysis_result, validation_errors = character(0),
           log = paste(log_lines, collapse = "\n"))
      
    }, error = function(e) {
      while (sink.number() > 0) sink()   # make sure output is restored
      log_lines <- if (file.exists(log_file)) readLines(log_file, warn = FALSE) else character(0)
      list(success = FALSE, error = conditionMessage(e), validation_errors = character(0),
           log = paste(log_lines, collapse = "\n"))
    })

    rv$run_output <- run_output

    if (isTRUE(run_output$success)) {
      res <- run_output$result
      rv$history <- rbind(rv$history, data.frame(
        Timestamp        = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        Target            = res$target_class_r,
        Comparison        = if (isTRUE(input$one_vs_all)) "ALL (one-vs-all)" else trimws(input$comparison_class),
        N_Runs            = round(input$n_runs),
        RF_AUC            = round(res$auc_df$Mean_AUC[1], 4),
        NN_PCA_AUC        = round(res$auc_df$Mean_AUC[2], 4),
        NN_ANOVA_AUC      = round(res$auc_df$Mean_AUC[3], 4),
        High_Confidence   = nrow(res$high_confidence_positives),
        Output_Folder     = output_dir_this_run,
        stringsAsFactors  = FALSE
      ))
    }

    # Jump the user to the Results tab automatically once a run finishes.
    if (isTRUE(run_output$success)) {
      updateNavbarPage(session, "main_nav", selected = "Results")
    }
  })

  # ---- Status line ----
  output$status <- renderText({
    ro <- rv$run_output
    if (isTRUE(rv$running)) return("Running - see progress bar...")
    if (is.null(ro)) return("Select your files (and check Settings), then click 'Run Analysis'.")
    if (!is.null(ro$validation_errors) && length(ro$validation_errors) > 0) {
      return("Please fix the issue(s) below before running.")
    }
    if (isTRUE(ro$success)) "Analysis complete - see the Results tab."
    else paste("ANALYSIS FAILED:", ro$error)
  })

  output$validation_errors <- renderUI({
    ro <- rv$run_output
    if (is.null(ro) || is.null(ro$validation_errors) || length(ro$validation_errors) == 0) {
      return(NULL)
    }
    tags$div(
      style = "background-color:#fdecea; border:1px solid #f5c2c0; padding:10px; border-radius:4px; margin-bottom:10px;",
      tags$strong("Please fix the following:"),
      tags$ul(lapply(ro$validation_errors, tags$li))
    )
  })

  has_results <- reactive({
    !is.null(rv$run_output) && isTRUE(rv$run_output$success)
  })

  # ---- Summary tab ----
  output$dataset_summary <- renderTable({
    req(has_results())
    rv$run_output$result$test_diagnostics
  })

  output$results_overview <- renderText({
    req(has_results())
    res <- rv$run_output$result
    paste0(
      "Target class              : ", res$target_class_r, "\n",
      "Descriptors selected      : ", length(res$descriptor_cols), "\n",
      "Total new peptides scored : ", nrow(res$simple_predictions), "\n",
      "High-confidence candidates: ", nrow(res$high_confidence_positives), "\n",
      "All 3 models agree (>=90%): ", sum(res$clear_predictions$Models_Agree_90 == 3)
    )
  })

  # ---- Model Performance tab ----
  output$model_performance <- renderTable({
    req(has_results())
    rv$run_output$result$auc_df
  })
  output$stability_table <- renderDT({
    req(has_results())
    datatable(rv$run_output$result$stability_summary, options = list(pageLength = 10))
  })

  # ---- Predictions tab ----
  output$simple_predictions_table <- renderDT({
    req(has_results())
    datatable(rv$run_output$result$simple_predictions, options = list(pageLength = 10))
  })
  output$detailed_predictions_table <- renderDT({
    req(has_results())
    datatable(rv$run_output$result$clear_predictions, options = list(pageLength = 10, scrollX = TRUE))
  })
  output$high_confidence_table <- renderDT({
    req(has_results())
    datatable(rv$run_output$result$high_confidence_positives, options = list(pageLength = 10, scrollX = TRUE))
  })
  output$consensus_table <- renderDT({
    req(has_results())
    datatable(rv$run_output$result$consensus, options = list(pageLength = 10, scrollX = TRUE))
  })

  # ---- Features tab ----
  output$importance_table <- renderDT({
    req(has_results())
    datatable(rv$run_output$result$importance_df, options = list(pageLength = 15))
  })
  output$feature_stats_table <- renderTable({
    req(has_results())
    rv$run_output$result$feature_stats
  })

  # ---- PCA tab ----
  # Served via Shiny's own resource path (addResourcePath) rather than
  # embedded as raw HTML via iframe srcdoc - the self-contained plotly
  # widgets are large and can fail to execute their scripts when
  # inlined as an srcdoc attribute in some browsers. Serving the actual
  # file over http:// is the same as how it behaves if you open it
  # directly, so it renders reliably.
  output$pca_train_view <- renderUI({
    req(has_results())
    pca_dir <- dirname(rv$run_output$result$pca_html_train_path)
    addResourcePath("pca_output_current", pca_dir)
    tags$iframe(
      src = file.path("pca_output_current", basename(rv$run_output$result$pca_html_train_path)),
      style = "width:100%; height:600px; border:none;"
    )
  })
  output$pca_combined_view <- renderUI({
    req(has_results())
    pca_dir <- dirname(rv$run_output$result$pca_html_combined_path)
    addResourcePath("pca_output_current", pca_dir)
    tags$iframe(
      src = file.path("pca_output_current", basename(rv$run_output$result$pca_html_combined_path)),
      style = "width:100%; height:600px; border:none;"
    )
  })

  # ---- Files / Output tab ----
  output$output_dir_text <- renderText({
    req(has_results())
    paste("Files were written to:",
          normalizePath(dirname(rv$run_output$result$excel_path), mustWork = FALSE))
  })

  output$download_buttons <- renderUI({
    req(has_results())
    tagList(
      downloadButton("download_excel", "Download Excel Workbook"),
      br(), br(),
      downloadButton("download_pca_train", "Download PCA (Training)"),
      br(), br(),
      downloadButton("download_pca_combined", "Download PCA (Training + New)")
    )
  })

  output$download_excel <- downloadHandler(
    filename = function() basename(rv$run_output$result$excel_path),
    content  = function(file) file.copy(rv$run_output$result$excel_path, file)
  )
  output$download_pca_train <- downloadHandler(
    filename = function() basename(rv$run_output$result$pca_html_train_path),
    content  = function(file) file.copy(rv$run_output$result$pca_html_train_path, file)
  )
  output$download_pca_combined <- downloadHandler(
    filename = function() basename(rv$run_output$result$pca_html_combined_path),
    content  = function(file) file.copy(rv$run_output$result$pca_html_combined_path, file)
  )

  # ---- Open Output Folder (desktop-style convenience) ----
  output$open_folder_button <- renderUI({
    req(has_results())
    actionButton("open_folder_btn", "Open Output Folder", icon = icon("folder-open"), width = "100%")
  })

  observeEvent(input$open_folder_btn, {
    req(has_results())
    folder <- normalizePath(dirname(rv$run_output$result$excel_path), mustWork = FALSE)
    tryCatch({
      if (.Platform$OS.type == "windows") {
        shell.exec(folder)
      } else {
        utils::browseURL(folder)
      }
    }, error = function(e) {
      showNotification(paste("Could not open the folder automatically. It is located at:", folder),
                        type = "warning", duration = 8)
    })
  })

  output$run_history_table <- renderDT({
    datatable(rv$history, options = list(pageLength = 5, scrollX = TRUE))
  })

  # ---- Log tab ----
  output$run_log <- renderText({
    ro <- rv$run_output
    if (is.null(ro)) return("(log will appear here after you click Run Analysis)")
    if (is.null(ro$log) || nchar(ro$log) == 0) return("(no log captured)")
    ro$log
  })
}

shinyApp(ui = ui, server = server)
