#
# LDH Cytotoxicity Assay Analysis Shiny App
# Version: 1.0
# Date: 2025-09-14
# Author: S. Berg
#
# This application processes, analyzes, and visualizes data from a standard
# LDH cytotoxicity plate reader assay. It supports per-treatment normalization
# and generates an interactive plot, statistical summaries, and a downloadable report.
#

# --- 1. Load Required Libraries ---
library(shiny)      # The core Shiny framework
library(tidyverse)  # For data manipulation (dplyr, tidyr) and plotting (ggplot2)
library(car)        # For Levene's test (assumption check)
library(broom)      # For tidying model outputs
library(ggpubr)     # For adding statistical annotations to plots
library(dunn.test)  # For Dunn's post-hoc test (non-parametric)
library(DT)         # For interactive data tables
library(ggiraph)    # For creating interactive ggplots
library(rmarkdown)  # For generating the downloadable report


# --- 2. Define the User Interface (UI) ---
ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "yeti"),

  titlePanel("LDH Cytotoxicity Assay Analysis"),

  sidebarLayout(
    # --- 2a. Sidebar Panel for Inputs ---
    sidebarPanel(
      width = 3,
      h4("1. Upload Data"),
      fileInput("layout_file", "Plate Layout File (CSV)", accept = ".csv"),
      fileInput("raw_data_file", "Raw Absorbance File (CSV)", accept = ".csv"),
      hr(),

      h4("2. Define Controls"),
      helpText("Select groups for calculation and comparison."),
      uiOutput("spontaneous_control_ui"),
      textInput("max_lysis_prefix", "Max Lysis Prefix", value = "Tritonx100-"),
      uiOutput("blank_control_ui"),
      hr(),

      h4("3. Analysis Settings"),
      uiOutput("comparison_control_ui"),
      radioButtons("stat_method_choice", "Statistical Test Method:",
                   choices = c("Automatic Selection" = "auto",
                               "Parametric (ANOVA)" = "anova",
                               "Non-parametric (Kruskal-Wallis)" = "kruskal"),
                   selected = "auto"),
      hr(),

      h4("4. Execute and Export"),
      actionButton("run_analysis", "Run Analysis", icon = icon("play-circle"), class = "btn-primary btn-lg"),
      downloadButton("download_report", "Download Report", class = "btn-success")
    ),

    # --- 2b. Main Panel for Outputs ---
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "results_tabs",
        tabPanel("Summary & Plot",
                 h3("Cytotoxicity Results", style = "margin-top:20px;"),
                 girafeOutput("cytotoxicity_plot", height = "500px"),
                 hr(),
                 h4("Analysis Summary"),
                 uiOutput("summary_text")),
        tabPanel("Detailed Statistics",
                 h3("Statistical Test Details", style = "margin-top:20px;"),
                 h4("Overall Test Result (Omnibus Test)"),
                 uiOutput("omnibus_test_results"),
                 hr(),
                 h4("Pairwise Comparisons (Post-Hoc Test)"),
                 DT::dataTableOutput("pairwise_table")),
        tabPanel("Processed Data Preview",
                 h3("Data after Normalization", style = "margin-top:20px"),
                 helpText("This table shows the calculated cytotoxicity data used for analysis."),
                 downloadButton("download_data", "Download Data as CSV", class = "btn-info"),
                 hr(),
                 DT::dataTableOutput("processed_data_table")),
        tabPanel("Assumption Checks",
                 uiOutput("assumption_checks_ui")),
        tabPanel("Interpretation Guide",
                 uiOutput("interpretation_guide_ui"))
      )
    )
  )
)


# --- 3. Define the Server Logic ---
server <- function(input, output) {

  # --- 3a. Dynamic UI Rendering ---
  # These reactive elements populate the dropdowns in the UI based on the uploaded layout file.

  layout_groups <- reactive({
    req(input$layout_file)
    df <- read.csv(input$layout_file$datapath, header = TRUE)
    df %>%
      pivot_longer(cols = everything(), names_to = "column", values_to = "group") %>%
      pull(group) %>%
      unique() %>%
      na.omit()
  })

  output$spontaneous_control_ui <- renderUI({
    req(layout_groups())
    selectInput("spontaneous_control", "Spontaneous Lysis (Untreated)",
                choices = layout_groups(), selected = "untreated")
  })

  output$blank_control_ui <- renderUI({
    req(layout_groups())
    selectInput("blank_control", "Blank (Media Only)",
                choices = layout_groups(), selected = "empty")
  })

  output$comparison_control_ui <- renderUI({
    req(layout_groups(), input$spontaneous_control, input$blank_control, input$max_lysis_prefix)
    all_groups <- layout_groups()
    calc_controls <- c(input$spontaneous_control, input$blank_control)
    triton_groups <- all_groups[str_detect(all_groups, fixed(input$max_lysis_prefix))]
    treatment_choices <- setdiff(all_groups, c(calc_controls, triton_groups))
    selectInput("comparison_control", "Comparison Control (e.g., Vehicle)", choices = treatment_choices)
  })


  # --- 3b. Main Analysis Reactive Block ---
  # This eventReactive block contains the entire analysis pipeline.
  # It only runs when the "Run Analysis" button is clicked.
  analysis_results <- eventReactive(input$run_analysis, {

    # 1. ===== INPUT VALIDATION =====
    # Ensure all necessary inputs are available before proceeding.
    req(input$layout_file, input$raw_data_file, input$spontaneous_control,
        input$max_lysis_prefix, input$comparison_control)

    plate_layout_raw <- try(read.csv(input$layout_file$datapath, header = TRUE))
    raw_absorbance <- try(read.csv(input$raw_data_file$datapath, header = FALSE))

    # Provide user-friendly error messages for common input problems.
    validate(
      need(inherits(plate_layout_raw, "data.frame"), "Error: Could not read Plate Layout file. Please ensure it is a valid CSV."),
      need(inherits(raw_absorbance, "data.frame"), "Error: Could not read Raw Absorbance file. Please ensure it is a valid CSV.")
    )

    plate_layout <- plate_layout_raw %>%
      mutate(Row = LETTERS[1:nrow(.)]) %>%
      pivot_longer(cols = -Row, names_to = "Column_Num", values_to = "group") %>%
      mutate(Column = as.numeric(gsub("X", "", Column_Num)), well = paste0(Row, Column)) %>%
      select(well, group)

    validate(
      need(input$spontaneous_control %in% plate_layout$group, "Error: The selected 'Spontaneous Lysis' group was not found."),
      need(input$blank_control %in% plate_layout$group, "Error: The selected 'Blank' group was not found."),
      need(any(str_detect(plate_layout$group, fixed(input$max_lysis_prefix))),
           paste0("Error: The prefix '", input$max_lysis_prefix, "' did not match any group names."))
    )

    # 2. ===== DATA PROCESSING =====
    # Background correction (680nm) and tidying data from wide to long format.
    abs_490 <- raw_absorbance %>% filter(row_number() %% 2 != 0)
    abs_680 <- raw_absorbance %>% filter(row_number() %% 2 == 0)
    corrected_abs <- abs_490 - abs_680

    corrected_long <- corrected_abs %>%
      mutate(Row = LETTERS[1:nrow(.)]) %>%
      pivot_longer(cols = -Row, names_to = "Column_Num", values_to = "absorbance") %>%
      mutate(Column = as.numeric(gsub("V", "", Column_Num)), well = paste0(Row, Column)) %>%
      select(well, absorbance)

    df_merged <- inner_join(plate_layout, corrected_long, by = "well")

    # 3. ===== NORMALIZATION CALCULATION =====
    # Calculate plate-wide blank and spontaneous lysis averages.
    blank_avg <- df_merged %>% filter(group == input$blank_control) %>% summarise(mean_abs = mean(absorbance, na.rm = TRUE)) %>% pull(mean_abs)
    df_blank_corrected <- df_merged %>% mutate(abs_corrected = absorbance - blank_avg)
    spontaneous_avg <- df_blank_corrected %>% filter(group == input$spontaneous_control) %>% summarise(mean_val = mean(abs_corrected, na.rm = TRUE)) %>% pull(mean_val)

    # Grouped normalization: calculate a unique max lysis for each treatment "family".
    df_with_family <- df_blank_corrected %>% mutate(treatment_family = str_remove(group, fixed(input$max_lysis_prefix)))
    df_max_lysis <- df_with_family %>% filter(str_detect(group, fixed(input$max_lysis_prefix))) %>% group_by(treatment_family) %>% summarise(maximum_avg = mean(abs_corrected, na.rm = TRUE))
    df_ready_for_calc <- df_with_family %>% left_join(df_max_lysis, by = "treatment_family")

    # 4. ===== FINAL CYTOTOXICITY CALCULATION =====
    # Apply the cytotoxicity formula and filter out all control groups.
    df_long <- df_ready_for_calc %>%
      mutate(
        denominator = maximum_avg - spontaneous_avg,
        cytotoxicity = ((abs_corrected - spontaneous_avg) / denominator) * 100,
        cytotoxicity = ifelse(cytotoxicity < 0, 0, cytotoxicity)
      ) %>%
      filter(
        !str_detect(group, fixed(input$max_lysis_prefix)),
        !group %in% c(input$spontaneous_control, input$blank_control),
        !str_detect(group, "ldh")
      ) %>%
      select(group, cytotoxicity) %>%
      mutate(group = factor(group)) # Re-factor to drop unused levels

    validate(
      need(nrow(df_long) > 0, "Error: No data remains after filtering controls."),
      need(all(is.finite(df_long$cytotoxicity)), "Calculation Error: Max Lysis may be <= Spontaneous Lysis.")
    )

    # 5. ===== STATISTICAL ANALYSIS =====
    control_group_name <- input$comparison_control
    model <- lm(cytotoxicity ~ group, data = df_long)

    # Assumption checks
    shapiro_result <- shapiro.test(residuals(model))
    levene_result <- leveneTest(cytotoxicity ~ group, data = df_long)

    # Determine statistical path (auto or manual)
    use_anova <- (shapiro_result$p.value > 0.05 && levene_result$`Pr(>F)`[1] > 0.05)
    if (input$stat_method_choice == "anova") use_anova <- TRUE
    if (input$stat_method_choice == "kruskal") use_anova <- FALSE

    # Execute the appropriate statistical test (parametric or non-parametric)
    if (use_anova) {
      test_name <- "One-Way ANOVA with Tukey's HSD"; stat_method_plot <- "t.test"; method_plot <- "anova"
      aov_model <- aov(model); anova_summary <- summary(aov_model); p_value <- anova_summary[[1]]$`Pr(>F)`[1]
      omnibus_text <- HTML(sprintf("The overall ANOVA test was <b>%s</b> (p-value = %.4f)", ifelse(p_value <= 0.05, "significant", "not significant"), p_value))
      pairwise_df <- as.data.frame(TukeyHSD(aov_model)$group) %>% tibble::rownames_to_column(var = "comparison") %>% rename(p.adjusted = `p adj`)
    } else {
      test_name <- "Kruskal-Wallis with Dunn's Test"; stat_method_plot <- "wilcox.test"; method_plot <- "kruskal.test"
      kruskal_result <- kruskal.test(cytotoxicity ~ group, data = df_long); p_value <- kruskal_result$p.value
      omnibus_text <- HTML(sprintf("The overall Kruskal-Wallis test was <b>%s</b> (p-value = %.4f)", ifelse(p_value <= 0.05, "significant", "not significant"), p_value))
      dunn_result <- dunn.test(df_long$cytotoxicity, df_long$group, method = "bonferroni")
      pairwise_df <- data.frame(comparison = dunn_result$comparisons, p.adjusted = dunn_result$P.adjusted)
    }

    # 6. ===== PLOT GENERATION =====
    # Create a list of pairs for ggpubr to draw significance brackets
    treatment_names <- as.character(unique(df_long$group))
    comparison_groups <- treatment_names[!treatment_names %in% control_group_name]
    comparisons_list <- lapply(comparison_groups, function(x) c(x, control_group_name))

    # Build the ggplot object with interactive layers
    final_plot <- ggplot(df_long, aes(x = group, y = cytotoxicity, fill = group)) +
      geom_boxplot_interactive(alpha = 0.7, outlier.shape = NA) +
      geom_jitter_interactive(
        width = 0.15, height = 0, size = 2.5,
        aes(color = group, tooltip = sprintf("Group: %s\nCytotoxicity: %.1f%%", group, cytotoxicity), data_id = group)
      ) +
      stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
      scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.05, 0.15))) +
      labs(title = "LDH Assay: Compound Cytotoxicity", x = "Treatment Group", y = "Cytotoxicity (%)") +
      theme_light(base_size = 14) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
      stat_compare_means(method = method_plot, label.y.npc = 0.95) +
      stat_compare_means(comparisons = comparisons_list, label = "p.signif", method = stat_method_plot)

    # 7. ===== RETURN RESULTS =====
    # Package all results into a list to be used by the output renderers.
    list(
      final_plot = final_plot,
      qq_plot = ggplot(broom::augment(model), aes(sample = .resid)) + geom_qq() + geom_qq_line() + theme_light(),
      assumption_check_text = paste("Shapiro-Wilk P-value:", round(shapiro_result$p.value, 4), "\nLevene's Test P-value:", round(levene_result$`Pr(>F)`[1], 4)),
      omnibus_text = omnibus_text,
      pairwise_df = pairwise_df,
      summary_text = HTML(paste0("<p>Analysis performed via <b>", test_name, "</b>.</p>")),
      processed_data = df_long,
      test_name = test_name,
      test_choice = input$stat_method_choice,
      omnibus_p_value = p_value,
      control_group_name = control_group_name,
      shapiro_p = shapiro_result$p.value,
      levene_p = levene_result$`Pr(>F)`[1]
    )
  })


  # --- 3c. Standard Output Rendering ---
  # These blocks render the primary outputs like plots and tables.
  output$cytotoxicity_plot <- renderGirafe({ req(analysis_results()); girafe(ggobj = analysis_results()$final_plot) })
  output$qq_plot <- renderPlot({ req(analysis_results()); analysis_results()$qq_plot })
  output$summary_text <- renderUI({ req(analysis_results()); analysis_results()$summary_text })
  output$omnibus_test_results <- renderUI({ req(analysis_results()); analysis_results()$omnibus_text })
  output$assumption_check_text <- renderPrint({ req(analysis_results()); cat(analysis_results()$assumption_check_text) })
  output$pairwise_table <- DT::renderDataTable({ req(analysis_results()); DT::datatable(analysis_results()$pairwise_df, options = list(pageLength = 10), rownames = FALSE) })
  output$processed_data_table <- DT::renderDataTable({ req(analysis_results()); DT::datatable(analysis_results()$processed_data, options = list(pageLength = 10), rownames = FALSE) })


  # --- 3d. Dynamic Tab Rendering ---
  # These blocks generate the dynamic UI for the interpretation and assumption check tabs.
  output$assumption_checks_ui <- renderUI({
    req(analysis_results())
    results <- analysis_results()
    normality_met <- results$shapiro_p > 0.05
    normality_text <- sprintf("The Shapiro-Wilk test checks for normality. Your data's p-value is <b>%.4f</b>. Since this is %s 0.05, the assumption of normality <b>%s</b>.",
                              results$shapiro_p, ifelse(normality_met, "greater than", "less than or equal to"), ifelse(normality_met, "is met", "is violated"))
    variance_met <- results$levene_p > 0.05
    variance_text <- sprintf("Levene's test checks for equal variances. Your data's p-value is <b>%.4f</b>. Since this is %s 0.05, the assumption of equal variances <b>%s</b>.",
                             results$levene_p, ifelse(variance_met, "greater than", "less than or equal to"), ifelse(variance_met, "is met", "is violated"))
    conclusion_text <- if (normality_met && variance_met) {"Because both assumptions were met, a <b>parametric test (ANOVA)</b> is the most appropriate choice."} else {"Because at least one assumption was violated, a <b>non-parametric test (Kruskal-Wallis)</b> is the safer choice."}

    tagList(
      h3("Diagnostic Checks for Parametric Tests", style = "margin-top:20px;"),
      p("Parametric tests like ANOVA are most reliable when data meets two key assumptions: normality and equal variances."),
      hr(), h4("1. Normality of Residuals"), HTML(normality_text),
      plotOutput("qq_plot", height = "400px"),
      helpText("The Q-Q plot is a visual check for normality. Points should fall along the line."),
      hr(), h4("2. Homogeneity of Variances"), HTML(variance_text),
      hr(), h4("Conclusion"), HTML(paste0("<p>", conclusion_text, "</p>")),
      hr(), h4("Numerical Test Results"), verbatimTextOutput("assumption_check_text")
    )
  })

  output$interpretation_guide_ui <- renderUI({
    req(analysis_results())
    results <- analysis_results()
    test_choice_text <- switch(results$test_choice, "auto" = "Based on assumption checks, the app automatically selected the", "anova" = "You manually selected the", "kruskal" = "You manually selected the")
    omnibus_sig_text <- ifelse(results$omnibus_p_value <= 0.05, "significant", "not significant")
    omnibus_conclusion <- ifelse(results$omnibus_p_value <= 0.05, "This suggests a difference exists somewhere among your groups.", "This suggests no significant differences were found.")
    
    example_text <- ""
    if (!is.null(results$control_group_name) && results$control_group_name != "") {
      example_row <- results$pairwise_df %>% filter(str_detect(comparison, fixed(results$control_group_name))) %>% slice(1)
      if (nrow(example_row) > 0) {
        comp_names <- example_row$comparison; comp_p_val <- example_row$p.adjusted; comp_sig <- ifelse(comp_p_val <= 0.05, "is", "is not"); comp_effect <- ifelse(comp_p_val <= 0.05, "a significant", "not a significant")
        example_text <- sprintf("<p>For example, the comparison between <b>%s</b> had a p-value of <b>%.4f</b>. As this value %s ≤ 0.05, we conclude there is %s difference between these groups.</p>",
                                comp_names, comp_p_val, comp_sig, comp_effect)
      }
    }

    tagList(
      h3("How to Interpret Your Results", style = "margin-top:20px;"),
      h4("1. The Overall (Omnibus) Test"),
      HTML(sprintf("<p>%s <b>%s</b>. This test gives one p-value to answer: 'Is there a difference anywhere among my groups?'</p><p>Your overall p-value was <b>%.4f</b>, which is <b>%s</b>.</p><p>%s</p>",
                   test_choice_text, results$test_name, results$omnibus_p_value, omnibus_sig_text, omnibus_conclusion)),
      hr(), h4("2. The Pairwise (Post-Hoc) Test"),
      HTML(paste0("<p>A post-hoc test was then run to compare every pair of groups. See 'Detailed Statistics' for the full table.</p>", example_text))
    )
  })


  # --- 3e. Download Handlers ---
  # These blocks handle the logic for the download buttons.
  output$download_data <- downloadHandler(
    filename = function() { paste("processed-cytotoxicity-data-", Sys.Date(), ".csv", sep = "") },
    content = function(file) {
      req(analysis_results())
      write.csv(analysis_results()$processed_data, file, row.names = FALSE)
    }
  )

  output$download_report <- downloadHandler(
    filename = function() { paste("LDH-Assay-Report-", Sys.Date(), ".html", sep = "") },
    content = function(file) {
      req(analysis_results())
      withProgress(message = 'Generating Report...', value = 0, {
        incProgress(0.2, detail = "Setting up...")
        tempReport <- file.path(tempdir(), "report.Rmd")
        file.copy("report.Rmd", tempReport, overwrite = TRUE) # Assumes report.Rmd is in the same directory as app.R
        report_params <- list(results = analysis_results())
        incProgress(0.6, detail = "Knitting document...")
        rmarkdown::render(tempReport, output_file = file, params = report_params, envir = new.env(parent = globalenv()))
      })
    }
  )

}


# --- 4. Run the Application ---
shinyApp(ui = ui, server = server)