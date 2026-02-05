# Packages (guard installs)
pkgs <- c("lmerTest","lme4","ggplot2","lattice","dplyr")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

# Data - set working directory to input file location
input_file  <- '/Users/allisonpickle/Desktop/Microglia 5ugLPS Mixed.csv'
setwd(dirname(input_file))
process_data <- read.csv(input_file)

cat("\nUnique Treatment values in raw CSV:\n")
print(sort(unique(process_data$Treatment)))

# Factors (hard-coded order, reused for plots)
baseline_treatment <- "Ctrl"   # baseline used for p-values and contrasts

# Create display labels for plotting with "-" and "+"
# The names of this vector define the specific order of treatment levels
treatment_labels <- c(
  "Ctrl"  = "-\n-",
  "5ug/mL"  = "-\n+\n"
)

# Extract levels directly from the names of the labels vector
treatment_levels <- names(treatment_labels)
process_data$Treatment <- factor(process_data$Treatment, levels = treatment_levels)

cat("\nTreatment counts after factor():\n")
print(table(process_data$Treatment, useNA = "ifany"))

# Define row labels for left side of strip plots
# Order corresponds to the rows in the strip plot (Top to Bottom)
row_labels <- c(
  "Control",
  "5ug/mL"
)

process_data$Well <- factor(process_data$Well)

# CHECK IF BIOLOGICAL_REP COLUMN EXISTS
if ("Biological_Rep" %in% names(process_data)) {
  process_data$Biological_Rep <- factor(process_data$Biological_Rep)
  message("Biological_Rep column found and converted to factor")
} else {
  stop("ERROR: No 'Biological_Rep' column found in CSV. Please check column names above and either:\n",
       "  1. Add a 'Biological_Rep' column to your CSV, or\n",
       "  2. Uncomment one of the options above to create/rename it")
}

# CHECK IF IMAGE COLUMN EXISTS (required for nested model)
if ("Image" %in% names(process_data)) {
  process_data$Image <- factor(process_data$Image)
  message("Image column found and converted to factor")
} else {
  warning("No 'Image' column found. Nested model will skip (1|Well:Image) term.")
}

# ----------------- MODELS (NESTED ONLY) -----------------

# Nested random effects model: (1|Biological_Rep) + (1|Biological_Rep:Well) + (1|Well:Image)
if ("Image" %in% names(process_data)) {
  NestedFit <- lmerTest::lmer(Results ~ Treatment + 
                                (1 | Biological_Rep) + 
                                (1 | Biological_Rep:Well) + 
                                (1 | Well:Image),
                              data = process_data)
  anova(NestedFit); summary(NestedFit)
} else {
  NestedFit <- lmerTest::lmer(Results ~ Treatment + 
                                (1 | Biological_Rep) + 
                                (1 | Biological_Rep:Well),
                              data = process_data)
  anova(NestedFit); summary(NestedFit)
}

# Null model for nested structure (LRT)
if ("Image" %in% names(process_data)) {
  NestedFit_null <- lme4::lmer(Results ~ 1 + 
                                 (1 | Biological_Rep) + 
                                 (1 | Biological_Rep:Well) + 
                                 (1 | Well:Image),
                               data = process_data,
                               REML = FALSE)
  NestedFit_full <- lme4::lmer(Results ~ Treatment + 
                                 (1 | Biological_Rep) + 
                                 (1 | Biological_Rep:Well) + 
                                 (1 | Well:Image),
                               data = process_data,
                               REML = FALSE)
} else {
  NestedFit_null <- lme4::lmer(Results ~ 1 + 
                                 (1 | Biological_Rep) + 
                                 (1 | Biological_Rep:Well),
                               data = process_data,
                               REML = FALSE)
  NestedFit_full <- lme4::lmer(Results ~ Treatment + 
                                 (1 | Biological_Rep) + 
                                 (1 | Biological_Rep:Well),
                               data = process_data,
                               REML = FALSE)
}

# Transform: log
process_data$logResults <- log(process_data$Results + 0.5)

# Nested model with log transformation
nested_formula <- formula(NestedFit)
log_nested_formula <- update(nested_formula, logResults ~ .)

LogNestedFit <- lmerTest::lmer(log_nested_formula, data = process_data)
anova(LogNestedFit); summary(LogNestedFit)

# Null model for log-transformed nested structure
if ("Image" %in% names(process_data)) {
  LogNestedFit_null <- lme4::lmer(logResults ~ 1 + 
                                    (1 | Biological_Rep) + 
                                    (1 | Biological_Rep:Well) + 
                                    (1 | Well:Image),
                                  data = process_data,
                                  REML = FALSE)
  LogNestedFit_full <- lme4::lmer(logResults ~ Treatment + 
                                    (1 | Biological_Rep) + 
                                    (1 | Biological_Rep:Well) + 
                                    (1 | Well:Image),
                                  data = process_data,
                                  REML = FALSE)
} else {
  LogNestedFit_null <- lme4::lmer(logResults ~ 1 + 
                                    (1 | Biological_Rep) + 
                                    (1 | Biological_Rep:Well),
                                  data = process_data,
                                  REML = FALSE)
  LogNestedFit_full <- lme4::lmer(logResults ~ Treatment + 
                                    (1 | Biological_Rep) + 
                                    (1 | Biological_Rep:Well),
                                  data = process_data,
                                  REML = FALSE)
}

# ----------------- ICC CALCULATION -----------------

compute_icc <- function(fit) {
  vc <- as.data.frame(VarCorr(fit))
  # Sum all random effect variances
  var_random <- sum(vc$vcov[vc$grp != "Residual"])
  # Residual variance
  var_resid  <- vc$vcov[vc$grp == "Residual"]
  # ICC = variance of random effects / total variance
  icc <- var_random / (var_random + var_resid)
  return(icc)
}

icc_nested <- compute_icc(NestedFit)
icc_log_nested <- compute_icc(LogNestedFit)

# ----------------- LIKELIHOOD RATIO TEST -----------------

lrt_nested <- anova(NestedFit_null, NestedFit_full)
lrt_log_nested <- anova(LogNestedFit_null, LogNestedFit_full)

# ----------------- DIAGNOSTIC PLOTS (interactive) -----------------

plot_diag <- function(fit, tag) {
  par(mfrow = c(1, 2))
  fitted_vals <- fitted(fit)
  std_resid   <- resid(fit, type = "pearson")
  plot(fitted_vals, std_resid,
       xlab = "Fitted Values",
       ylab = "Standardized (Pearson) Residuals",
       main = paste(tag, "Std Residuals vs Fitted"))
  abline(h = 0, col = "red", lwd = 2)
  
  qq_data <- qqnorm(resid(fit), plot.it = FALSE)
  plot(qq_data$x, qq_data$y, 
       xlab = "Theoretical Quantiles", 
       ylab = "Sample Quantiles",
       main = paste(tag, "QQ Plot"),
       pch = 16, col = "black")
  qqline(resid(fit), col = "red", lwd = 2)
  par(mfrow = c(1, 1))
}

# ----------------- SAVE MODEL SUMMARIES -----------------

save_model_summaries <- function() {
  csv_base <- tools::file_path_sans_ext(basename(input_file))
  out_dir  <- paste0("MM_Outputs_", csv_base)
  dir.create(out_dir, showWarnings = FALSE)

  out_file <- file.path(out_dir, "Mixed_Model_Summaries.txt")

  sink(out_file)
  cat("================================================================\n")
  cat("MIXED EFFECTS MODEL SUMMARIES (NESTED MODELS ONLY)\n")
  cat("================================================================\n\n")
  
  nested_formula <- formula(NestedFit)
  log_nested_formula <- formula(LogNestedFit)
  
  cat("=== NestedFit Model (Raw Data) ===\n")
  cat("Formula: "); print(nested_formula)
  cat("\nThis model accounts for:\n")
  cat("  - Treatment effects (fixed)\n")
  cat("  - Nested random effects:\n")
  cat("    * Biological replicate variation: (1|Biological_Rep)\n")
  cat("    * Well within biological replicate: (1|Biological_Rep:Well)\n")
  if ("Image" %in% names(process_data)) {
    cat("    * Image within well: (1|Well:Image)\n")
  }
  cat("\n")
  
  cat("ANOVA:\n"); print(anova(NestedFit))
  cat("\nSummary:\n"); print(summary(NestedFit))
  
  cat("\n\nRandom Effects Variance Components:\n")
  cat("----------------------------------------------------------------\n")
  vc_nested <- as.data.frame(VarCorr(NestedFit))
  print(vc_nested)
  
  cat("\n\nIntraclass Correlation Coefficient (ICC):\n")
  cat("----------------------------------------------------------------\n")
  cat(sprintf("ICC = %.4f\n", icc_nested))
  cat("Interpretation: Proportion of total variance attributable to all random effects.\n")
  
  cat("\n\nLikelihood Ratio Test (LRT) - Treatment Effect:\n")
  cat("----------------------------------------------------------------\n")
  print(lrt_nested)

  cat("\n\n================================================================\n")
  cat("=== LogNestedFit Model (Log-Transformed Data) ===\n")
  cat("Formula: "); print(log_nested_formula)
  cat("\nThis model accounts for:\n")
  cat("  - Treatment effects (fixed)\n")
  cat("  - Nested random effects with log-transformed response\n")
  cat("  - Transform applied: logResults = log(Results + 0.5)\n\n")
  
  cat("ANOVA:\n"); print(anova(LogNestedFit))
  cat("\nSummary:\n"); print(summary(LogNestedFit))
  
  cat("\n\nRandom Effects Variance Components:\n")
  cat("----------------------------------------------------------------\n")
  vc_log_nested <- as.data.frame(VarCorr(LogNestedFit))
  print(vc_log_nested)
  
  cat("\n\nIntraclass Correlation Coefficient (ICC):\n")
  cat("----------------------------------------------------------------\n")
  cat(sprintf("ICC = %.4f\n", icc_log_nested))
  cat("Interpretation: Proportion of total variance attributable to all random effects.\n")
  
  cat("\n\nLikelihood Ratio Test (LRT) - Treatment Effect:\n")
  cat("----------------------------------------------------------------\n")
  print(lrt_log_nested)
  
  cat("\n\n================================================================\n")
  cat("Model Structure Notes:\n")
  cat("----------------------------------------------------------------\n")
  cat("NestedFit uses:    "); cat(deparse(nested_formula), "\n")
  cat("LogNestedFit uses: "); cat(deparse(log_nested_formula), "\n")
  cat("\nBoth models include hierarchical random effects structure:\n")
  cat("  Biological_Rep > Biological_Rep:Well")
  if ("Image" %in% names(process_data)) {
    cat(" > Well:Image")
  }
  cat("\n\nThis structure properly accounts for:\n")
  cat("  1. Between-biological-replicate variation\n")
  cat("  2. Between-well variation within the same biological replicate\n")
  if ("Image" %in% names(process_data)) {
    cat("  3. Between-image variation within the same well\n")
  }
  cat("\n")
  
  sink()

  message("Saved ANOVA and summary outputs to: ", out_file)
}

# ----------------- SAVE DIAGNOSTIC IMAGES & RESIDUALS -----------------

save_diagnostics <- function(save_data = TRUE) {
  csv_base <- tools::file_path_sans_ext(basename(input_file))
  out_dir  <- paste0("MM_Outputs_", csv_base)
  dir.create(out_dir, showWarnings = FALSE)
  
  models <- list(
    NestedFit = NestedFit,
    LogNestedFit = LogNestedFit
  )
  
  all_resid <- list()
  
  for (nm in names(models)) {
    fit <- models[[nm]]
    fitted_vals <- fitted(fit)
    std_resid   <- resid(fit, type = "pearson")
    
    png(file.path(out_dir, paste0(nm, "_residuals.png")),
        width = 1200, height = 1000, res = 150)
    plot(fitted_vals, std_resid,
         xlab = "Fitted Values",
         ylab = "Standardized (Pearson) Residuals",
         main = paste(nm, "Std Residuals vs Fitted"))
    abline(h = 0, col = "red", lwd = 2)
    dev.off()
    
    png(file.path(out_dir, paste0(nm, "_qq.png")),
        width = 1200, height = 1000, res = 150)
    qq_data <- qqnorm(resid(fit), plot.it = FALSE)
    plot(qq_data$x, qq_data$y, 
         xlab = "Theoretical Quantiles", 
         ylab = "Sample Quantiles",
         main = paste(nm, "QQ Plot"),
         pch = 16, col = "#0068f9")
    qqline(resid(fit), col = "red", lwd = 2)
    dev.off()
    
    if (save_data) {
      all_resid[[nm]] <- data.frame(
        Model       = nm,
        Fitted      = fitted_vals,
        StdResidual = std_resid,
        RawResidual = resid(fit)
      )
    }
  }
  
  if (save_data) {
    write.csv(
      do.call(rbind, all_resid),
      file.path(out_dir, "All_Calculated_Residuals.csv"),
      row.names = FALSE
    )
  }
  
  message("Saved standardized (Pearson) residual vs fitted and QQ plots for each model to: ", out_dir)
}

# ----------------- P-VALUES AND BARPLOTS -----------------

extract_treatment_pvals <- function(fit, baseline_name = "Ctrl") {
  fe <- coef(summary(fit))
  rn <- rownames(fe)
  is_trt <- grepl("^Treatment", rn)
  if (!any(is_trt)) {
    stop("No Treatment* rows found in fixed effects.")
  }

  trt_rows  <- fe[is_trt, , drop = FALSE]
  trt_names <- sub("^Treatment", "", rownames(trt_rows))

  p_col <- grep("Pr\\(", colnames(trt_rows), value = TRUE)[1]
  if (is.na(p_col)) stop("Could not find p-value column in fixed effects.")

  pvals <- trt_rows[, p_col]
  names(pvals) <- trt_names

  data.frame(
    Treatment = c(baseline_name, trt_names),
    p.value   = c(NA, as.numeric(pvals)),
    stringsAsFactors = FALSE
  )
}

add_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

build_bar_summary_with_sig <- function(data, fit, response_col, baseline_name = baseline_treatment) {
  summary_df <- data |>
    dplyr::group_by(Treatment) |>
    dplyr::summarise(
      mean = mean(.data[[response_col]], na.rm = TRUE),
      se   = sd(.data[[response_col]], na.rm = TRUE) /
             sqrt(sum(!is.na(.data[[response_col]]))),
      .groups = "drop"
    )

  summary_df$Treatment <- factor(
    summary_df$Treatment,
    levels = levels(data$Treatment)
  )

  p_df <- extract_treatment_pvals(fit, baseline_name)
  p_df$Treatment <- factor(p_df$Treatment, levels = levels(data$Treatment))

  out <- dplyr::left_join(summary_df, p_df, by = "Treatment")
  out$stars <- vapply(out$p.value, add_stars, character(1))
  out
}

plot_bar_p_values <- function(summ_df, raw_data, response_col, y_label, title_text) {
  max_y <- max(summ_df$mean + summ_df$se, na.rm = TRUE)
  
  aligned_treatment_labels <- setNames(paste0("\n", treatment_labels), names(treatment_labels))
  combined_row_label <- paste0("\n", paste(row_labels, collapse = "\n"))
  
  ggplot(summ_df, aes(x = Treatment, y = mean)) +
    geom_col(fill = "#6c6c6c", alpha = 0.7) +
    geom_jitter(data = raw_data, 
                aes(x = Treatment, y = .data[[response_col]], shape = Biological_Rep, color = Biological_Rep), 
                width = 0.2, 
                alpha = 0.6,
                size = 1.0) +
    scale_shape_manual(values = c(16, 17, 15, 18, 3, 4, 8, 7)) +
    scale_color_manual(values = c("#08267d", "#228B22", "#FF8C00", "black", "black", "black", "black", "black")) +
    scale_x_discrete(labels = aligned_treatment_labels) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    geom_text(
      aes(y = mean + se, label = stars),
      vjust = -0.5,
      size = 7,
      fontface = "bold"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
    labs(
      x = "",
      y = y_label, 
      title = title_text,
      shape = "Biological Rep",
      color = "Biological Rep"
    ) +
    annotate("text", x = 0.4, y = -Inf, label = combined_row_label, 
             hjust = 1, vjust = 1.1, size = 3, color = "black", lineheight = 0.8) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, color = "black", lineheight = 0.8, size = 9, margin = margin(t = 0)),
      axis.text.y = element_text(color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, color = "black"),
      plot.margin = margin(5, 5, 15, 50),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    ) +
    coord_cartesian(clip = "off")
}

# ----------------- BUILD AND SAVE BARPLOTS -----------------

csv_base <- tools::file_path_sans_ext(basename(input_file))
out_dir  <- paste0("MM_Outputs_", csv_base)
dir.create(out_dir, showWarnings = FALSE)

# NestedFit barplot (raw data)
nested_bar_summary <- build_bar_summary_with_sig(
  data         = process_data,
  fit          = NestedFit,
  response_col = "Results",
  baseline_name = baseline_treatment
)

nested_bar_plot <- plot_bar_p_values(
  nested_bar_summary,
  raw_data     = process_data,
  response_col = "Results",
  y_label      = "Results (mean ± SE)",
  title_text   = "Treatment (NestedFit)"
)

ggplot2::ggsave(
  filename = file.path(out_dir, "Barplot_NestedFit.png"),
  plot     = nested_bar_plot,
  width    = 6,
  height   = 4,
  dpi      = 300
)

# LogNestedFit barplot (log-transformed data)
log_nested_bar_summary <- build_bar_summary_with_sig(
  data         = process_data,
  fit          = LogNestedFit,
  response_col = "logResults",
  baseline_name = baseline_treatment
)

log_nested_bar_plot <- plot_bar_p_values(
  log_nested_bar_summary,
  raw_data     = process_data,
  response_col = "logResults",
  y_label      = "log(Results + 0.5) (mean ± SE)",
  title_text   = "Treatment (LogNestedFit)"
)

ggplot2::ggsave(
  filename = file.path(out_dir, "Barplot_LogNestedFit.png"),
  plot     = log_nested_bar_plot,
  width    = 6,
  height   = 4,
  dpi      = 300
)

# ---- Run saves at the end ----
save_diagnostics()
save_model_summaries()