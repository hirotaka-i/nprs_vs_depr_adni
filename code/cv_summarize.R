#!/usr/bin/env Rscript
# cv_summarize.R
# Summarize cross-validation results and recommend parameters
# Usage: Rscript cv_summarize.R <phenotype.pheno>

library(tidyverse)
library(pROC)
library(meta)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  cat("Usage: Rscript cv_summarize.R <phenotype.pheno>\n")
  quit(status = 1)
}

pheno_file <- args[1]
outcome <- pheno_file %>% basename() %>% str_replace("\\.pheno$", "")
cv_dir <- paste0("temp/cv5/", outcome)
cat(cv_dir)
results_dir <- file.path(cv_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

cat("========================================\n")
cat("Summarizing CV Results\n")
cat("========================================\n\n")

# ============================================================================
# Load data
# ============================================================================

cat("Loading training results...\n")
training_files <- list.files(file.path(cv_dir, "training"), 
                              pattern = "fold_.*\\.prsice$", 
                              full.names = TRUE)

train_results <- map_df(training_files, function(f) {
  basename <- basename(f)
  fold <- as.integer(str_extract(basename, "(?<=fold_)\\d+"))
  r2_str <- str_extract(basename, "(?<=r2_)[0-9_]+")
  r2 <- as.numeric(gsub("_", ".", r2_str))
  
  data <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
  data$fold <- fold
  data$clump_r2 <- r2
  data$phase <- "train"
  return(data)
})

cat("Loading test results...\n")
test_results_file <- file.path(cv_dir, "testing", "test_results.csv")
test_results <- read.csv(test_results_file, stringsAsFactors = FALSE)

best_params <- read.csv(file.path(cv_dir, "training", "best_params_per_fold.csv"),
                        stringsAsFactors = FALSE)

cat(sprintf("Loaded %d training results across %d folds\n", 
            nrow(train_results), length(unique(train_results$fold))))
cat(sprintf("Loaded %d test results\n\n", nrow(test_results)))

# ============================================================================
# Test performance summary
# ============================================================================

cat("========================================\n")
cat("Test Performance Summary\n")
cat("========================================\n\n")

test_summary <- test_results %>%
  summarise(
    n_folds = n(),
    mean_R2 = mean(test_R2),
    se_R2 = sd(test_R2) / sqrt(n()),
    ci_lower_R2 = mean_R2 - 1.96 * se_R2,
    ci_upper_R2 = mean_R2 + 1.96 * se_R2,
    min_R2 = min(test_R2),
    max_R2 = max(test_R2)
  )

cat(sprintf("Mean test R² = %.4f ± %.4f (SE)\n", 
            test_summary$mean_R2, test_summary$se_R2))
cat(sprintf("95%% CI: [%.4f, %.4f]\n", 
            test_summary$ci_lower_R2, test_summary$ci_upper_R2))
cat(sprintf("Range: [%.4f, %.4f]\n\n", 
            test_summary$min_R2, test_summary$max_R2))

valid_p <- test_results$test_P[test_results$test_P > 0]
if (length(valid_p) > 0) {
  geom_mean_p <- exp(mean(log(valid_p)))
  cat(sprintf("Geometric mean P-value = %.2e\n", geom_mean_p))
  cat(sprintf("Median P-value = %.2e\n", median(valid_p)))
  cat(sprintf("Significant folds (P < 0.05): %d / %d\n\n", 
              sum(test_results$test_P < 0.05), nrow(test_results)))
}

# ============================================================================
# Aggregate results by parameter combination
# ============================================================================

cat("========================================\n")
cat("Performance by Parameter Combination\n")
cat("========================================\n\n")

train_aggregated <- train_results %>%
  group_by(clump_r2, Threshold) %>%
  summarise(
    mean_R2_train = mean(R2),
    se_R2_train = sd(R2) / sqrt(n()),
    mean_P_train = exp(mean(log(pmax(P, 1e-300)))),
    .groups = "drop"
  )

best_counts <- best_params %>%
  count(best_r2, best_p, name = "n_folds_selected")

param_summary <- test_results %>%
  group_by(test_r2, test_p) %>%
  summarise(
    mean_R2_test = mean(test_R2),
    se_R2_test = sd(test_R2) / sqrt(n()),
    mean_P_test = exp(mean(log(pmax(test_P, 1e-300)))),
    mean_num_snp = mean(num_snp),
    .groups = "drop"
  ) %>%
  rename(clump_r2 = test_r2, p_threshold = test_p) %>%
  left_join(best_counts, by = c("clump_r2" = "best_r2", "p_threshold" = "best_p")) %>%
  left_join(train_aggregated, by = c("clump_r2" = "clump_r2", "p_threshold" = "Threshold")) %>%
  replace_na(list(n_folds_selected = 0)) %>%
  arrange(desc(mean_R2_test))

print(param_summary)

write.csv(param_summary, 
          file.path(results_dir, "cv_results_aggregated.csv"),
          row.names = FALSE, quote = FALSE)

# ============================================================================
# Recommend parameters
# ============================================================================

cat("\n========================================\n")
cat("Parameter Recommendations\n")
cat("========================================\n\n")

modal_params <- best_params %>%
  count(best_r2, best_p, sort = TRUE) %>%
  slice(1)

cat("MODAL RECOMMENDATION (most frequently selected):\n")
cat(sprintf("  clump-r2 = %s\n", modal_params$best_r2))
cat(sprintf("  p-threshold = %s\n", modal_params$best_p))
cat(sprintf("  Selected in %d / 5 folds\n", modal_params$n))

modal_test <- test_results %>%
  filter(test_r2 == modal_params$best_r2, test_p == modal_params$best_p) %>%
  summarise(mean_R2 = mean(test_R2), se_R2 = sd(test_R2) / sqrt(n()))

if (nrow(modal_test) > 0) {
  cat(sprintf("  Test R² = %.4f ± %.4f\n", modal_test$mean_R2, modal_test$se_R2))
}

best_avg_params <- param_summary %>%
  filter(n_folds_selected > 0) %>%
  arrange(desc(mean_R2_test)) %>%
  slice(1)

cat("\nBEST AVERAGE RECOMMENDATION (highest mean test R²):\n")
cat(sprintf("  clump-r2 = %s\n", best_avg_params$clump_r2))
cat(sprintf("  p-threshold = %s\n", best_avg_params$p_threshold))
cat(sprintf("  Mean test R² = %.4f ± %.4f\n", 
            best_avg_params$mean_R2_test, best_avg_params$se_R2_test))
cat(sprintf("  Selected in %d / 5 folds\n\n", best_avg_params$n_folds_selected))

# Save outputs
write.csv(train_results, 
          file.path(results_dir, "cv_results_full.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(test_results, 
          file.path(results_dir, "cv_results_test.csv"),
          row.names = FALSE, quote = FALSE)

# ============================================================================
# Generate Visualization Figures
# ============================================================================

cat("========================================\n")
cat("Generating Visualization Figures\n")
cat("========================================\n\n")

pheno_data <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)

# ============================================================================
# Figure 1: Forest Plot
# ============================================================================

cat("Creating Figure 1: Forest plot...\n")

forest_data <- data.frame()

for (fold_id in 1:5) {
  prs_file <- file.path(cv_dir, "testing", paste0("fold_", fold_id, "_test.best"))
  
  if (!file.exists(prs_file)) {
    cat(sprintf("  Warning: PRS file not found for fold %d\n", fold_id))
    next
  }
  
  prs_scores <- read.table(prs_file, header = TRUE, stringsAsFactors = FALSE)
  
  fold_data <- pheno_data %>%
    inner_join(prs_scores %>% select(IID, PRS), by = "IID")
  
  fold_data$PRS_std <- scale(fold_data$PRS)[,1]
  
  formula_str <- "PHENO ~ PRS_std + AGE + PTGENDER + APOE4 + PTEDUCAT + DX_MCI + DX_AD + PC1 + PC2 + PC3 + PC4"
  model <- glm(as.formula(formula_str), data = fold_data, family = binomial())
  
  coef_summary <- summary(model)$coefficients["PRS_std", ]
  
  forest_data <- rbind(forest_data, data.frame(
    fold = fold_id,
    beta = coef_summary["Estimate"],
    se = coef_summary["Std. Error"],
    p = coef_summary["Pr(>|z|)"],
    OR = exp(coef_summary["Estimate"]),
    OR_lower = exp(coef_summary["Estimate"] - 1.96 * coef_summary["Std. Error"]),
    OR_upper = exp(coef_summary["Estimate"] + 1.96 * coef_summary["Std. Error"])
  ))
}

if (nrow(forest_data) > 0) {
  meta_result <- metagen(TE = beta, seTE = se, data = forest_data, 
                         sm = "OR", method.tau = "DL")
  
  forest_data <- rbind(forest_data, data.frame(
    fold = "Meta",
    beta = meta_result$TE.fixed,
    se = meta_result$seTE.fixed,
    p = meta_result$pval.fixed,
    OR = exp(meta_result$TE.fixed),
    OR_lower = exp(meta_result$lower.fixed),
    OR_upper = exp(meta_result$upper.fixed)
  ))
  
  # Create forest plot with proper margins
  png(file.path(results_dir, "figure1_forest_plot.png"), 
      width = 12, height = 6, units = "in", res = 300)
  
  par(mar = c(5, 6, 4, 8))
  
  y_pos <- nrow(forest_data):1
  xlim_range <- c(min(forest_data$OR_lower) * 0.8, max(forest_data$OR_upper) * 1.3)
  
  plot(forest_data$OR, y_pos, 
       xlim = xlim_range, ylim = c(0.5, nrow(forest_data) + 0.5),
       xlab = "Odds Ratio (95% CI) per SD increase in PRS", ylab = "",
       main = "Forest Plot: Neuroticism PRS and Depression Risk",
       pch = ifelse(forest_data$fold == "Meta", 18, 16),
       cex = ifelse(forest_data$fold == "Meta", 1.5, 1.2),
       col = ifelse(forest_data$fold == "Meta", "red", "black"),
       axes = FALSE)
  
  for (i in 1:nrow(forest_data)) {
    lines(c(forest_data$OR_lower[i], forest_data$OR_upper[i]), 
          c(y_pos[i], y_pos[i]), 
          lwd = ifelse(forest_data$fold[i] == "Meta", 2, 1),
          col = ifelse(forest_data$fold[i] == "Meta", "red", "black"))
  }
  
  abline(v = 1, lty = 2, col = "gray50")
  
  axis(1)
  axis(2, at = y_pos, labels = paste("Fold", forest_data$fold), las = 1)
  
  # Add text annotations outside plot
  text(par("usr")[2], y_pos, 
       sprintf("%.2f [%.2f, %.2f] P=%.3f", 
               forest_data$OR, forest_data$OR_lower, forest_data$OR_upper, 
               forest_data$p),
       pos = 4, cex = 0.8, xpd = TRUE)
  
  dev.off()
  
  cat("  Figure 1 saved: figure1_forest_plot.png\n")
  cat(sprintf("  Meta-analysis OR = %.3f (95%% CI: %.3f-%.3f), P = %.4f\n",
              exp(meta_result$TE.fixed), exp(meta_result$lower.fixed), 
              exp(meta_result$upper.fixed), meta_result$pval.fixed))
}

# ============================================================================
# Figure 2: Pooled ROC Curve
# ============================================================================

cat("\nCreating Figure 2: Pooled ROC curve...\n")

pooled_data <- data.frame()

for (fold_id in 1:5) {
  prs_file <- file.path(cv_dir, "testing", paste0("fold_", fold_id, "_test.best"))
  
  if (!file.exists(prs_file)) next
  
  prs_scores <- read.table(prs_file, header = TRUE, stringsAsFactors = FALSE)
  
  fold_data <- pheno_data %>%
    inner_join(prs_scores %>% select(IID, PRS), by = "IID") %>%
    mutate(fold = fold_id)
  
  fold_data$PRS_std <- scale(fold_data$PRS)[,1]
  
  pooled_data <- rbind(pooled_data, fold_data)
}

if (nrow(pooled_data) > 0) {
  roc_obj <- roc(pooled_data$PHENO, pooled_data$PRS_std, quiet = TRUE)
  
  png(file.path(results_dir, "figure2_roc_curve.png"), 
      width = 8, height = 8, units = "in", res = 300)
  
  plot(roc_obj, 
       main = sprintf("ROC Curve: Pooled Test Sets\nAUC = %.3f (95%% CI: %.3f-%.3f)",
                      auc(roc_obj), ci.auc(roc_obj)[1], ci.auc(roc_obj)[3]),
       col = "blue", lwd = 2,
       print.auc = FALSE)
  
  abline(a = 0, b = 1, lty = 2, col = "gray50")
  
  dev.off()
  
  cat("  Figure 2 saved: figure2_roc_curve.png\n")
  cat(sprintf("  Pooled AUC = %.3f (95%% CI: %.3f-%.3f)\n",
              auc(roc_obj), ci.auc(roc_obj)[1], ci.auc(roc_obj)[3]))
}

# ============================================================================
# Figure 3: ROC Curves by Fold with Mean and 95% CI
# ============================================================================
cat("\nCreating Figure 3: ROC curves by fold...\n")

# Initialize list for ROC objects
roc_list <- list()
auc_values <- numeric(5)

# Loop through each fold
for (fold_id in 1:5) {
  prs_file <- file.path(cv_dir, "testing", paste0("fold_", fold_id, "_test.best"))
  
  if (!file.exists(prs_file)) {
    cat(sprintf("  Warning: PRS file not found for fold %d\n", fold_id))
    next
  }
  
  prs_scores <- read.table(prs_file, header = TRUE, stringsAsFactors = FALSE)
  
  fold_data <- pheno_data %>%
    inner_join(prs_scores %>% select(IID, PRS), by = "IID")
  
  if (nrow(fold_data) == 0) {
    cat(sprintf("  Warning: No data for fold %d after merge\n", fold_id))
    next
  }
  
  fold_data$PHENO <- as.numeric(as.character(fold_data$PHENO))
  
  if (length(unique(fold_data$PHENO)) < 2) {
    cat(sprintf("  Warning: Fold %d lacks binary outcome variance\n", fold_id))
    next
  }
  
  fold_data$PRS_std <- scale(fold_data$PRS)[,1]
  
  roc_obj <- roc(response = fold_data$PHENO,
                 predictor = fold_data$PRS_std,
                 levels = c(0, 1),
                 quiet = TRUE)
  
  roc_list[[fold_id]] <- roc_obj
  auc_values[fold_id] <- auc(roc_obj)
  cat(sprintf("  Fold %d: AUC = %.3f\n", fold_id, auc(roc_obj)))
}

# Remove NULL entries
roc_list <- roc_list[!sapply(roc_list, is.null)]
auc_values <- auc_values[auc_values > 0]

if (length(roc_list) < 2) {
  cat("  ERROR: Fewer than 2 valid ROC curves. Skipping Figure 3.\n")
} else {
  # Calculate mean AUC with 95% CI
  mean_auc <- mean(auc_values)
  se_auc <- sd(auc_values) / sqrt(length(auc_values))
  ci_lower_auc <- mean_auc - 1.96 * se_auc
  ci_upper_auc <- mean_auc + 1.96 * se_auc
  
  # Manual calculation of pointwise mean and CI
  specificities_grid <- seq(0, 1, by = 0.01)
  
  sensitivities_list <- lapply(roc_list, function(roc_obj) {
    coords(roc_obj, x = specificities_grid, input = "specificity", 
           ret = "sensitivity")$sensitivity
  })
  
  sensitivities_matrix <- do.call(rbind, sensitivities_list)
  
  mean_sensitivities <- colMeans(sensitivities_matrix)
  ci_lower_sens <- apply(sensitivities_matrix, 2, quantile, probs = 0.025)
  ci_upper_sens <- apply(sensitivities_matrix, 2, quantile, probs = 0.975)
  
  # Create plot
  png(file.path(results_dir, "figure3_roc_by_fold.png"), 
      width = 10, height = 10, units = "in", res = 300)
  
  # Setup plot area
  plot(NULL, xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate (1 - Specificity)",
       ylab = "True Positive Rate (Sensitivity)",
       main = sprintf("ROC Curves: 5-Fold Cross-Validation\nMean AUC = %.3f (95%% CI: %.3f-%.3f)",
                      mean_auc, ci_lower_auc, ci_upper_auc),
       asp = 1)
  
  # Add individual ROC curves - MANUALLY EXTRACT COORDINATES
  for (roc_obj in roc_list) {
    # Extract sensitivities and specificities
    sens <- roc_obj$sensitivities
    spec <- roc_obj$specificities
    fpr <- 1 - spec  # Convert to FPR
    
    # Plot using FPR (not specificity)
    lines(fpr, sens, col = adjustcolor("black", alpha.f = 0.7), 
          lwd = 1, lty = "dashed")
  }
  
  # Add chance line
  abline(a = 0, b = 1, col = "red", lty = 2, lwd = 1.5)
  
  # Convert specificity grid to FPR
  fpr_grid <- 1 - rev(specificities_grid)
  
  # Plot 95% CI polygon
  polygon(x = c(fpr_grid, rev(fpr_grid)),
          y = c(rev(ci_lower_sens), ci_upper_sens),
          col = adjustcolor("steelblue", alpha.f = 0.3),
          border = NA)
  
  # Plot mean ROC line
  lines(x = fpr_grid, y = rev(mean_sensitivities), col = "darkblue", lwd = 3)
  
  # Add legend
  legend("bottomright",
         legend = c(sprintf("Mean ROC (AUC = %.3f)", mean_auc),
                    "95% Confidence Band",
                    "Individual Fold ROCs",
                    "Random Chance"),
         col = c("darkblue", adjustcolor("steelblue", alpha.f = 0.3), "gray60", "red"),
         lwd = c(3, 10, 1, 1.5),
         lty = c(1, NA, 2, 2),
         fill = c(NA, adjustcolor("steelblue", alpha.f = 0.3), NA, NA),
         border = NA,
         cex = 0.9,
         bty = "n")
  
  dev.off()
  
  cat("  Figure 3 saved: figure3_roc_by_fold.png\n")
  cat(sprintf("  Mean AUC = %.3f ± %.3f (SE)\n", mean_auc, se_auc))
}

# ============================================================================
# Figure 4: Quantile Analysis (3 panels)
# ============================================================================

cat("\nCreating Figure 4: Quantile analysis...\n")

quantile_data <- data.frame()

for (fold_id in 1:5) {
  prs_file <- file.path(cv_dir, "testing", paste0("fold_", fold_id, "_test.best"))
  
  if (!file.exists(prs_file)) next
  
  prs_scores <- read.table(prs_file, header = TRUE, stringsAsFactors = FALSE)
  
  fold_data <- pheno_data %>%
    inner_join(prs_scores %>% select(IID, PRS), by = "IID") %>%
    mutate(fold = fold_id)
  
  fold_data$PRS_std <- scale(fold_data$PRS)[,1]
  
  fold_data$quantile <- cut(fold_data$PRS_std, 
                             breaks = quantile(fold_data$PRS_std, probs = seq(0, 1, 0.25)),
                             labels = c("Q1", "Q2", "Q3", "Q4"),
                             include.lowest = TRUE)
  
  quantile_data <- rbind(quantile_data, 
                         fold_data %>% select(IID, PHENO, PRS_std, quantile, fold, 
                                              AGE, PTGENDER, APOE4, PTEDUCAT, 
                                              DX_MCI, DX_AD, PC1, PC2, PC3, PC4))
}

if (nrow(quantile_data) > 0) {
  # Calculate prevalence summary
  quantile_summary <- quantile_data %>%
    group_by(quantile) %>%
    summarise(
      n = n(),
      n_cases = sum(PHENO == 1),
      prevalence = mean(PHENO == 1),
      se = sqrt(prevalence * (1 - prevalence) / n),
      .groups = "drop"
    ) %>%
    mutate(
      ci_lower = pmax(0, prevalence - 1.96 * se),
      ci_upper = pmin(1, prevalence + 1.96 * se)
    )
  
  # Fit logistic regression model with quantile as predictor
  # Set Q1 as reference category
  quantile_data$quantile <- relevel(factor(quantile_data$quantile), ref = "Q1")
  
  model <- glm(PHENO ~ quantile + AGE + PTGENDER + APOE4 + PTEDUCAT + 
                 DX_MCI + DX_AD + PC1 + PC2 + PC3 + PC4,
               data = quantile_data,
               family = binomial())
  
  # Extract coefficients for quantiles
  coef_summary <- summary(model)$coefficients
  quantile_coefs <- coef_summary[grep("^quantile", rownames(coef_summary)), , drop = FALSE]
  
  # Create coefficient data frame (add Q1 as reference)
  coef_data <- data.frame(
    quantile = c("Q1", "Q2", "Q3", "Q4"),
    beta = c(0, quantile_coefs[, "Estimate"]),
    se = c(0, quantile_coefs[, "Std. Error"]),
    p = c(1, quantile_coefs[, "Pr(>|z|)"])
  ) %>%
    mutate(
      ci_lower = beta - 1.96 * se,
      ci_upper = beta + 1.96 * se,
      OR = exp(beta),
      OR_lower = exp(ci_lower),
      OR_upper = exp(ci_upper)
    )
  
  # Create 3-panel plot
  png(file.path(results_dir, "figure4_quantile_analysis.png"), 
      width = 10, height = 15, units = "in", res = 300)
  
  par(mfrow = c(3, 1), mar = c(4, 4, 3, 2))
  
  # Panel A: PRS distribution by quantile
  boxplot(PRS_std ~ quantile, data = quantile_data,
          main = "A. PRS Distribution by Quantile",
          xlab = "PRS Quantile", ylab = "Standardized PRS",
          col = c("#E8F4F8", "#B3D9E8", "#7EBED8", "#4AA3C8"),
          border = "darkblue")
  abline(h = 0, lty = 2, col = "gray50")
  
  # Panel B: Depression prevalence by quantile
  bp <- barplot(quantile_summary$prevalence, 
                names.arg = quantile_summary$quantile,
                main = "B. Depression Prevalence by PRS Quantile (Unadjusted)",
                xlab = "PRS Quantile", ylab = "Depression Prevalence",
                col = c("#E8F4F8", "#B3D9E8", "#7EBED8", "#4AA3C8"),
                border = "darkblue",
                ylim = c(0, max(quantile_summary$ci_upper) * 1.3))
  
  arrows(bp, quantile_summary$ci_lower, 
         bp, quantile_summary$ci_upper,
         angle = 90, code = 3, length = 0.1, lwd = 2)
  
  text(bp, max(quantile_summary$ci_upper) * 1.2,
       sprintf("n=%d\n(%d cases)", quantile_summary$n, quantile_summary$n_cases),
       cex = 0.8, pos = 1, xpd = TRUE)
  
  # Panel C: Adjusted log odds ratios (Beta coefficients)
  plot(1:4, coef_data$beta,
       xlim = c(0.5, 4.5), 
       ylim = range(c(coef_data$ci_lower, coef_data$ci_upper)) * 1.2,
       xlab = "PRS Quantile", ylab = "Log Odds Ratio (Beta)",
       main = "C. Adjusted Association: Depression Risk by PRS Quantile (vs Q1)",
       pch = 19, cex = 1.5, col = c("#E8F4F8", "#B3D9E8", "#7EBED8", "#4AA3C8"),
       xaxt = "n")
  
  # Add error bars
  arrows(1:4, coef_data$ci_lower, 
         1:4, coef_data$ci_upper,
         angle = 90, code = 3, length = 0.1, lwd = 2,
         col = c("#E8F4F8", "#B3D9E8", "#7EBED8", "#4AA3C8"))
  
  # Add reference line at 0
  abline(h = 0, lty = 2, col = "red")
  
  # Add x-axis labels
  axis(1, at = 1:4, labels = c("Q1\n(Ref)", "Q2", "Q3", "Q4"))
  
  # Add p-values
  text(1:4, coef_data$ci_upper + 0.05,
       sprintf("P=%.3f", coef_data$p),
       cex = 0.8, pos = 3)
  
  # Add OR text on right side
  text(rep(4.3, 4), coef_data$beta,
       sprintf("OR=%.2f\n[%.2f-%.2f]", 
               coef_data$OR, coef_data$OR_lower, coef_data$OR_upper),
       cex = 0.7, pos = 4, xpd = TRUE)
  
  dev.off()
  
  cat("  Figure 4 saved: figure4_quantile_analysis.png\n")
  cat("\n  Prevalence by quantile:\n")
  print(quantile_summary)
  cat("\n  Adjusted coefficients:\n")
  print(coef_data)
}

# Save supplementary data
write.csv(forest_data, 
          file.path(results_dir, "forest_plot_data.csv"),
          row.names = FALSE, quote = FALSE)

write.csv(data.frame(fold = 1:5, AUC = auc_values), 
          file.path(results_dir, "auc_by_fold.csv"),
          row.names = FALSE, quote = FALSE)

write.csv(quantile_summary, 
          file.path(results_dir, "quantile_summary.csv"),
          row.names = FALSE, quote = FALSE)

write.csv(coef_data,
          file.path(results_dir, "quantile_regression_coefs.csv"),
          row.names = FALSE, quote = FALSE)

cat("\n========================================\n")
cat("Summary Complete!\n")
cat("========================================\n")
cat("Files saved to:", results_dir, "\n")
cat("  - cv_summary.txt\n")
cat("  - recommended_params.txt\n")
cat("  - cv_results_aggregated.csv\n")
cat("  - cv_results_test.csv\n")
cat("  - cv_results_full.csv\n")
cat("\nFigures:\n")
cat("  - figure1_forest_plot.png\n")
cat("  - figure2_roc_curve.png\n")
cat("  - figure3_auc_by_fold.png\n")
cat("  - figure4_quantile_boxplot.png\n")