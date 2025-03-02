# Zero-Inflation Explorer: Interactive Dashboard for viromePower
# This Shiny app provides an interactive interface for exploring zero-inflation patterns
# and their impact on statistical power in virome studies.

library(shiny)
library(shinydashboard)
library(ggplot2)
library(plotly)
library(DT)
library(viromePower)

# UI Definition
ui <- dashboardPage(
  dashboardHeader(title = "Zero-Inflation Explorer"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview", tabName = "overview", icon = icon("dashboard")),
      menuItem("Simulation", tabName = "simulation", icon = icon("random")),
      menuItem("Parameter Estimation", tabName = "estimation", icon = icon("calculator")),
      menuItem("Power Analysis", tabName = "power", icon = icon("chart-line")),
      menuItem("Model Comparison", tabName = "comparison", icon = icon("balance-scale")),
      menuItem("Sampling Strategy", tabName = "sampling", icon = icon("search")),
      menuItem("Help", tabName = "help", icon = icon("question-circle"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # Overview Tab
      tabItem(tabName = "overview",
        fluidRow(
          box(
            title = "Welcome to Zero-Inflation Explorer",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            p("This interactive dashboard helps you explore zero-inflation patterns in virome data and analyze their impact on statistical power and study design."),
            p("Use the sidebar menu to navigate through different analysis modules:"),
            tags$ul(
              tags$li(tags$b("Simulation:"), "Generate and visualize simulated virome data with different zero-inflation patterns"),
              tags$li(tags$b("Parameter Estimation:"), "Estimate zero-inflation parameters from real virome data"),
              tags$li(tags$b("Power Analysis:"), "Evaluate statistical power for different zero-inflation scenarios"),
              tags$li(tags$b("Model Comparison:"), "Compare fixed vs. variable zero-inflation modeling approaches"),
              tags$li(tags$b("Sampling Strategy:"), "Optimize sample size and sequencing depth for given zero-inflation patterns")
            ),
            p("Start by exploring the simulation tab to understand how different zero-inflation patterns look in virome data.")
          )
        ),
        fluidRow(
          box(
            title = "About Variable Zero-Inflation",
            status = "info",
            solidHeader = TRUE,
            width = 6,
            p("Zero-inflation in virome data arises from two sources:"),
            tags$ul(
              tags$li(tags$b("Structural zeros:"), "True absence of a virus in a sample"),
              tags$li(tags$b("Sampling zeros:"), "Virus is present but not detected due to technical limitations")
            ),
            p("Traditional approaches use a single zero-inflation rate for all viral taxa, but in reality:"),
            tags$ul(
              tags$li("Different viruses have different prevalence patterns"),
              tags$li("Some viruses are genuinely rare while others are common but hard to detect"),
              tags$li("Zero-inflation rates may vary with abundance or other biological factors")
            ),
            p("Variable zero-inflation modeling uses a Beta distribution to capture this heterogeneity, resulting in more realistic simulations and more accurate power analysis.")
          ),
          box(
            title = "Beta Distribution Patterns",
            status = "info",
            solidHeader = TRUE,
            width = 6,
            plotOutput("betaPatterns")
          )
        )
      ),
      
      # Simulation Tab
      tabItem(tabName = "simulation",
        fluidRow(
          box(
            title = "Simulation Parameters",
            status = "primary",
            solidHeader = TRUE,
            width = 4,
            numericInput("n_samples", "Number of samples per group:", 30, min = 5, max = 200),
            numericInput("n_viruses", "Number of viral taxa:", 100, min = 10, max = 1000),
            numericInput("mean_zi", "Mean zero-inflation rate:", 0.7, min = 0.1, max = 0.9, step = 0.05),
            numericInput("sampling_zeros", "Sampling zero rate:", 0.2, min = 0, max = 0.5, step = 0.05),
            numericInput("effect_size", "Effect size (fold change):", 2.0, min = 1.1, max = 5.0, step = 0.1),
            checkboxInput("variable_zi", "Use variable zero-inflation", TRUE),
            conditionalPanel(
              condition = "input.variable_zi == true",
              numericInput("zi_alpha", "Beta shape parameter alpha:", 2, min = 0.1, max = 10, step = 0.1),
              numericInput("zi_beta", "Beta shape parameter beta:", 5, min = 0.1, max = 10, step = 0.1)
            ),
            actionButton("simulate", "Generate Simulation", class = "btn-primary")
          ),
          tabBox(
            title = "Simulation Results",
            width = 8,
            tabPanel("Zero-Inflation Distribution", 
                     plotlyOutput("ziDistPlot"),
                     downloadButton("downloadZIDist", "Download Plot")),
            tabPanel("Abundance vs. ZI", 
                     plotlyOutput("ziAbundancePlot"),
                     downloadButton("downloadZIAbundance", "Download Plot")),
            tabPanel("Group Comparison", 
                     plotlyOutput("ziGroupPlot"),
                     downloadButton("downloadZIGroup", "Download Plot")),
            tabPanel("Count Matrix", 
                     DT::dataTableOutput("countTable"),
                     downloadButton("downloadCounts", "Download Data"))
          )
        )
      ),
      
      # Parameter Estimation Tab
      tabItem(tabName = "estimation",
        fluidRow(
          box(
            title = "Upload Virome Data",
            status = "primary",
            solidHeader = TRUE,
            width = 4,
            fileInput("dataFile", "Upload count matrix (CSV):",
                      accept = c("text/csv", "text/comma-separated-values", ".csv")),
            numericInput("presenceThreshold", "Presence threshold:", 1, min = 0, max = 100),
            numericInput("minPrevalence", "Minimum prevalence:", 0.05, min = 0, max = 1, step = 0.01),
            selectInput("estimationMethod", "Estimation method:",
                       choices = c("moment", "mle", "bayes"),
                       selected = "mle"),
            actionButton("estimate", "Estimate Parameters", class = "btn-primary")
          ),
          tabBox(
            title = "Estimation Results",
            width = 8,
            tabPanel("Fitted Distribution", 
                     plotlyOutput("fittedDistPlot"),
                     verbatimTextOutput("fittedParams"),
                     downloadButton("downloadFitPlot", "Download Plot")),
            tabPanel("Diagnostic Plots", 
                     plotlyOutput("diagnosticPlot"),
                     downloadButton("downloadDiagPlot", "Download Plot")),
            tabPanel("Parameter Details", 
                     verbatimTextOutput("paramDetails"),
                     downloadButton("downloadParams", "Download Parameters"))
          )
        )
      ),
      
      # Power Analysis Tab
      tabItem(tabName = "power",
        fluidRow(
          box(
            title = "Power Analysis Settings",
            status = "primary",
            solidHeader = TRUE,
            width = 4,
            numericInput("power_n_samples", "Sample size per group:", 30, min = 5, max = 200),
            numericInput("power_effect_size", "Effect size:", 2.0, min = 1.1, max = 5.0, step = 0.1),
            numericInput("power_n_viruses", "Number of viral taxa:", 100, min = 10, max = 1000),
            numericInput("power_mean_zi", "Mean zero-inflation rate:", 0.7, min = 0.1, max = 0.9, step = 0.05),
            numericInput("power_sampling_zeros", "Sampling zero rate:", 0.2, min = 0, max = 0.5, step = 0.05),
            selectInput("power_zi_patterns", "Zero-inflation patterns to compare:",
                       choices = c("fixed", "right_skewed", "left_skewed", "symmetric", "uniform"),
                       selected = c("fixed", "right_skewed"),
                       multiple = TRUE),
            numericInput("power_n_sim", "Number of simulations:", 5, min = 5, max = 100, step = 5),
            actionButton("runPower", "Run Power Analysis", class = "btn-primary")
          ),
          tabBox(
            title = "Power Analysis Results",
            width = 8,
            tabPanel("Power Comparison", 
                     plotlyOutput("powerComparisonPlot"),
                     downloadButton("downloadPowerComp", "Download Plot")),
            tabPanel("Sample Size Effects", 
                     plotlyOutput("sampleSizePlot"),
                     downloadButton("downloadSampleSize", "Download Plot")),
            tabPanel("Effect Size Impact", 
                     plotlyOutput("effectSizePlot"),
                     downloadButton("downloadEffectSize", "Download Plot")),
            tabPanel("Summary Table", 
                     DT::dataTableOutput("powerSummaryTable"),
                     downloadButton("downloadPowerSummary", "Download Summary"))
          )
        )
      ),
      
      # Model Comparison Tab
      tabItem(tabName = "comparison",
        fluidRow(
          box(
            title = "Model Comparison Settings",
            status = "primary",
            solidHeader = TRUE,
            width = 4,
            fileInput("modelDataFile", "Upload count matrix (CSV):",
                      accept = c("text/csv", "text/comma-separated-values", ".csv")),
            numericInput("model_threshold", "Presence threshold:", 1, min = 0, max = 100),
            selectInput("model_method", "Variable ZI method:",
                       choices = c("beta", "power_law", "mixed"),
                       selected = "beta"),
            numericInput("model_bootstrap", "Bootstrap samples:", 100, min = 10, max = 1000),
            checkboxInput("model_visualize", "Generate visualizations", TRUE),
            actionButton("compareModels", "Compare Models", class = "btn-primary")
          ),
          tabBox(
            title = "Model Comparison Results",
            width = 8,
            tabPanel("Model Fit", 
                     plotlyOutput("modelFitPlot"),
                     verbatimTextOutput("modelFitSummary"),
                     downloadButton("downloadModelFit", "Download Plot")),
            tabPanel("ZI Distribution", 
                     plotlyOutput("modelZIPlot"),
                     downloadButton("downloadModelZI", "Download Plot")),
            tabPanel("Abundance Relationship", 
                     plotlyOutput("modelAbundancePlot"),
                     downloadButton("downloadModelAbundance", "Download Plot")),
            tabPanel("Statistical Tests", 
                     verbatimTextOutput("modelTestResults"),
                     downloadButton("downloadModelTests", "Download Results"))
          )
        )
      ),
      
      # Sampling Strategy Tab
      tabItem(tabName = "sampling",
        fluidRow(
          box(
            title = "Sampling Strategy Optimization",
            status = "primary",
            solidHeader = TRUE,
            width = 4,
            numericInput("strategy_target_power", "Target power:", 0.8, min = 0.5, max = 0.95, step = 0.05),
            numericInput("strategy_effect_size", "Effect size:", 2.0, min = 1.1, max = 5.0, step = 0.1),
            numericInput("strategy_mean_zi", "Mean zero-inflation rate:", 0.7, min = 0.1, max = 0.9, step = 0.05),
            numericInput("strategy_zi_alpha", "Beta shape parameter alpha:", 2, min = 0.1, max = 10, step = 0.1),
            numericInput("strategy_zi_beta", "Beta shape parameter beta:", 5, min = 0.1, max = 10, step = 0.1),
            numericInput("strategy_sequencing_effort", "Total sequencing effort (million reads):", 100, min = 10, max = 1000),
            numericInput("strategy_depth_min", "Min sequencing depth (million reads):", 0.5, min = 0.1, max = 10),
            numericInput("strategy_depth_max", "Max sequencing depth (million reads):", 5, min = 1, max = 20),
            checkboxInput("strategy_test_variable", "Compare fixed vs. variable ZI", TRUE),
            actionButton("optimizeStrategy", "Optimize Strategy", class = "btn-primary")
          ),
          tabBox(
            title = "Optimization Results",
            width = 8,
            tabPanel("Recommended Strategy", 
                     verbatimTextOutput("strategyRecommendation"),
                     downloadButton("downloadRecommendation", "Download")),
            tabPanel("Power Contour", 
                     plotlyOutput("strategyContourPlot"),
                     downloadButton("downloadContour", "Download Plot")),
            tabPanel("Comparison Plot", 
                     plotlyOutput("strategyComparisonPlot"),
                     downloadButton("downloadStrategyComp", "Download Plot")),
            tabPanel("Efficiency Plot", 
                     plotlyOutput("strategyEfficiencyPlot"),
                     downloadButton("downloadEfficiency", "Download Plot"))
          )
        )
      ),
      
      # Help Tab
      tabItem(tabName = "help",
        fluidRow(
          box(
            title = "Using the Zero-Inflation Explorer",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            p("This dashboard helps researchers understand and account for zero-inflation patterns in virome data analysis."),
            tags$h4("Quick Start Guide"),
            tags$ol(
              tags$li("Begin in the", tags$b("Simulation"), "tab to explore how different zero-inflation patterns look"),
              tags$li("Use the", tags$b("Parameter Estimation"), "tab to infer zero-inflation patterns from your own data"),
              tags$li("Analyze statistical power with the", tags$b("Power Analysis"), "tab"),
              tags$li("Compare fixed vs. variable zero-inflation models in the", tags$b("Model Comparison"), "tab"),
              tags$li("Optimize your study design using the", tags$b("Sampling Strategy"), "tab")
            ),
            tags$h4("Understanding Zero-Inflation Parameters"),
            tags$ul(
              tags$li(tags$b("Mean zero-inflation rate:"), "The average proportion of structural zeros across all viral taxa"),
              tags$li(tags$b("Sampling zero rate:"), "The proportion of detection failures among non-structural zeros"),
              tags$li(tags$b("Alpha and Beta:"), "Shape parameters for the Beta distribution that controls the distribution of zero-inflation rates across viral taxa:"),
              tags$ul(
                tags$li(tags$i("Alpha < Beta (e.g., 2, 5):"), "Right-skewed distribution - more viruses with lower ZI rates"),
                tags$li(tags$i("Alpha > Beta (e.g., 5, 2):"), "Left-skewed distribution - more viruses with higher ZI rates"),
                tags$li(tags$i("Alpha = Beta (e.g., 2, 2):"), "Symmetric distribution centered at 0.5")
              )
            ),
            tags$h4("Downloading Results"),
            p("Each tab includes download buttons for saving plots and data for use in presentations or publications."),
            tags$h4("Additional Help"),
            p("For more information on the viromePower package and its functions, refer to the package documentation:"),
            tags$code("?viromePower"),
            br(), br(),
            p("If you encounter any issues or have questions, please submit them to:"),
            tags$a(href="https://github.com/username/viromePower/issues", "https://github.com/username/viromePower/issues")
          )
        )
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  
  # Reactive values to store simulation and analysis results
  sim_data <- reactiveVal(NULL)
  param_estimates <- reactiveVal(NULL)
  power_results <- reactiveVal(NULL)
  model_comparison <- reactiveVal(NULL)
  sampling_strategy <- reactiveVal(NULL)
  
  # Overview tab - Show Beta distribution patterns
  output$betaPatterns <- renderPlot({
    x <- seq(0, 1, by = 0.01)
    df <- data.frame(
      x = rep(x, 4),
      y = c(
        dbeta(x, 2, 5),  # Right-skewed
        dbeta(x, 5, 2),  # Left-skewed
        dbeta(x, 2, 2),  # Symmetric
        dbeta(x, 1, 1)   # Uniform
      ),
      Pattern = rep(c("Right-skewed (2, 5)", "Left-skewed (5, 2)", 
                    "Symmetric (2, 2)", "Uniform (1, 1)"), each = length(x))
    )
    
    ggplot(df, aes(x = x, y = y, color = Pattern)) +
      geom_line(size = 1.2) +
      labs(
        x = "Zero-Inflation Rate",
        y = "Probability Density",
        title = "Beta Distribution Patterns for Zero-Inflation"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  # Simulation tab - Generate data when button is clicked
  observeEvent(input$simulate, {
    # Show progress
    withProgress(message = "Generating simulation...", {
      # Generate simulation based on input parameters
      sim <- simulate_zero_inflated_virome(
        n_samples = input$n_samples * 2,  # Total samples across both groups
        n_viruses = input$n_viruses,
        structural_zeros = input$mean_zi,
        sampling_zeros = input$sampling_zeros,
        effect_size = input$effect_size,
        variable_zi_rates = input$variable_zi,
        zi_alpha = if(input$variable_zi) input$zi_alpha else NA,
        zi_beta = if(input$variable_zi) input$zi_beta else NA
      )
      
      # Store the simulation
      sim_data(sim)
    })
  })
  
  # Render Zero-Inflation Distribution Plot
  output$ziDistPlot <- renderPlotly({
    req(sim_data())
    
    # Use the plot_virus_specific_zi_rates function and convert to plotly
    p <- plot_virus_specific_zi_rates(sim_data())
    ggplotly(p)
  })
  
  # Render Abundance vs. ZI Plot
  output$ziAbundancePlot <- renderPlotly({
    req(sim_data())
    
    # Use the plot_zi_by_abundance function and convert to plotly
    p <- plot_zi_by_abundance(sim_data())
    ggplotly(p)
  })
  
  # Render Group Comparison Plot
  output$ziGroupPlot <- renderPlotly({
    req(sim_data())
    
    # Use the compare_group_zi_patterns function and convert to plotly
    p <- compare_group_zi_patterns(sim_data(), plot_type = "scatter")
    ggplotly(p)
  })
  
  # Render Count Matrix Table
  output$countTable <- DT::renderDataTable({
    req(sim_data())
    
    # Format count data for display
    counts <- sim_data()$counts
    
    # Return as datatable
    DT::datatable(
      counts,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip'
      ),
      rownames = TRUE,
      caption = "Simulated Count Matrix (rows = viral taxa, columns = samples)"
    )
  })
  
  # Parameter Estimation tab - Estimate parameters when button is clicked
  observeEvent(input$estimate, {
    req(input$dataFile)
    
    # Show progress
    withProgress(message = "Estimating parameters...", {
      # Read uploaded data
      count_matrix <- read.csv(input$dataFile$datapath, row.names = 1)
      
      # Estimate parameters
      estimates <- estimate_zi_parameters_from_data(
        count_matrix = as.matrix(count_matrix),
        presence_threshold = input$presenceThreshold,
        min_prevalence = input$minPrevalence,
        method = input$estimationMethod,
        plot_fit = TRUE
      )
      
      # Store the estimates
      param_estimates(estimates)
    })
  })
  
  # Render Fitted Distribution Plot
  output$fittedDistPlot <- renderPlotly({
    req(param_estimates())
    
    # Get the plot from estimation results
    p <- param_estimates()$plot
    ggplotly(p)
  })
  
  # Render Parameter Details
  output$fittedParams <- renderPrint({
    req(param_estimates())
    
    # Extract key parameters
    params <- param_estimates()
    cat("Estimated Zero-Inflation Parameters:\n")
    cat("--------------------------------\n")
    cat("Alpha:", round(params$alpha, 3), "\n")
    cat("Beta:", round(params$beta, 3), "\n")
    cat("Mean ZI rate:", round(params$mean_zi, 3), "\n")
    cat("Median ZI rate:", round(params$median_zi, 3), "\n")
    cat("SD of ZI rates:", round(params$sd_zi, 3), "\n")
    cat("Range of ZI rates:", round(params$min_zi, 3), "to", round(params$max_zi, 3), "\n")
  })
  
  # Power Analysis tab - Run power analysis when button is clicked
  observeEvent(input$runPower, {
    # Show progress
    withProgress(message = "Running power analysis...", {
      # Create zi_patterns based on selections
      selected_patterns <- input$power_zi_patterns
      zi_patterns <- list()
      
      if ("fixed" %in% selected_patterns) {
        zi_patterns$fixed <- c(alpha = NA, beta = NA)
      }
      if ("right_skewed" %in% selected_patterns) {
        zi_patterns$right_skewed <- c(alpha = 2, beta = 5)
      }
      if ("left_skewed" %in% selected_patterns) {
        zi_patterns$left_skewed <- c(alpha = 5, beta = 2)
      }
      if ("symmetric" %in% selected_patterns) {
        zi_patterns$symmetric <- c(alpha = 2, beta = 2)
      }
      if ("uniform" %in% selected_patterns) {
        zi_patterns$uniform <- c(alpha = 1, beta = 1)
      }
      
      # Run power analysis across zi patterns
      power_result <- plot_power_by_zi_pattern(
        sample_sizes = c(input$power_n_samples * 0.5, input$power_n_samples, input$power_n_samples * 1.5),
        effect_sizes = c(input$power_effect_size * 0.75, input$power_effect_size, input$power_effect_size * 1.25),
        n_viruses = input$power_n_viruses,
        zi_patterns = zi_patterns,
        n_sim = input$power_n_sim,
        fixed_zi_rate = input$power_mean_zi,
        sampling_zeros = input$power_sampling_zeros
      )
      
      # Store the results
      power_results(power_result)
    })
  })
  
  # Render Power Comparison Plot
  output$powerComparisonPlot <- renderPlotly({
    req(power_results())
    
    # Get combined plot from power analysis results
    p <- power_results()$plots$combined
    # Convert to plotly if it's a ggplot
    if (inherits(p, "ggplot")) {
      ggplotly(p)
    } else {
      # For grid.arrange objects, we need to render as is
      output$powerComparisonPlot <- renderPlot({
        p
      })
      NULL
    }
  })
  
  # Model Comparison tab - Run model comparison when button is clicked
  observeEvent(input$compareModels, {
    req(input$modelDataFile)
    
    # Show progress
    withProgress(message = "Comparing models...", {
      # Read uploaded data
      count_matrix <- read.csv(input$modelDataFile$datapath, row.names = 1)
      
      # Run model comparison
      comparison <- compare_zi_models(
        count_matrix = as.matrix(count_matrix),
        presence_threshold = input$model_threshold,
        variable_zi_method = input$model_method,
        n_bootstrap = input$model_bootstrap,
        visualize = input$model_visualize
      )
      
      # Store the results
      model_comparison(comparison)
    })
  })
  
  # Render Model Fit Plot
  output$modelFitPlot <- renderPlotly({
    req(model_comparison())
    
    # Get model comparison plot
    p <- model_comparison()$plots$model_comparison
    if (is.null(p)) {
      return(NULL)
    }
    ggplotly(p)
  })
  
  # Render Model Fit Summary
  output$modelFitSummary <- renderPrint({
    req(model_comparison())
    
    # Extract key results
    results <- model_comparison()
    cat("Model Comparison Results:\n")
    cat("--------------------------------\n")
    cat("Best model:", results$best_model, "\n")
    cat("Likelihood ratio:", round(results$likelihood_ratio, 3), "\n")
    cat("LRT p-value:", formatC(results$lrt_p_value, digits = 4, format = "f"), "\n")
    cat("\nAIC values:\n")
    cat("Fixed model:", round(results$aic_values$fixed, 2), "\n")
    cat("Variable model:", round(results$aic_values$variable, 2), "\n")
    cat("Difference (Fixed - Variable):", round(results$aic_values$difference, 2), "\n")
    
    cat("\nBIC values:\n")
    cat("Fixed model:", round(results$bic_values$fixed, 2), "\n")
    cat("Variable model:", round(results$bic_values$variable, 2), "\n")
    cat("Difference (Fixed - Variable):", round(results$bic_values$difference, 2), "\n")
  })
  
  # Sampling Strategy tab - Run optimization when button is clicked
  observeEvent(input$optimizeStrategy, {
    # Show progress
    withProgress(message = "Optimizing sampling strategy...", {
      # Run optimization
      strategy <- optimize_sampling_strategy(
        target_power = input$strategy_target_power,
        effect_size = input$strategy_effect_size,
        zi_alpha = input$strategy_zi_alpha,
        zi_beta = input$strategy_zi_beta,
        total_sequencing_effort = input$strategy_sequencing_effort,
        mean_zi = input$strategy_mean_zi,
        depth_range = c(input$strategy_depth_min, input$strategy_depth_max),
        test_variable_zi = input$strategy_test_variable
      )
      
      # Store the results
      sampling_strategy(strategy)
    })
  })
  
  # Render Strategy Recommendation
  output$strategyRecommendation <- renderPrint({
    req(sampling_strategy())
    
    # Output the recommendation text
    cat(sampling_strategy()$recommendation)
  })
  
  # Render Strategy Contour Plot
  output$strategyContourPlot <- renderPlotly({
    req(sampling_strategy())
    
    # Get appropriate contour plot
    if (input$strategy_test_variable) {
      p <- sampling_strategy()$plots$contour_variable
    } else {
      p <- sampling_strategy()$plots$contour_plot
    }
    
    # Convert to plotly if it's a ggplot
    if (inherits(p, "ggplot")) {
      ggplotly(p)
    } else {
      NULL
    }
  })
  
  # Download handlers for all plots and data
  # (implementations would go here)
  
}

# Run the application
shinyApp(ui = ui, server = server)