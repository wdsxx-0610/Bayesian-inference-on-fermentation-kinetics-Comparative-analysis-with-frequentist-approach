# ==============================================================================
# Visualization: Frequentist vs Bayesian Predictions
# ==============================================================================

# --- 1. Load packages ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, cowplot, scales, readr)

# --- 2. Set paths (matching previous run results) ---
# Prediction results folder
result_dir <- "./results/Half_prediction_Results"
# Original data file path
raw_data_path <- "../data/CofermentationData.csv"

# --- 3. Read original data ---
if(file.exists(raw_data_path)){
  real_data <- read.csv(raw_data_path) %>%
    rename(time = time) %>% # Ensure time column is named time
    rename_with(~ str_to_title(.), .cols = -time) # Ensure column names are capitalized (Lactate, etc.)
  
  # Fix PDO column name (if it became Pdo -> PDO)
  if ("Pdo" %in% colnames(real_data)) {
    real_data <- real_data %>% rename(PDO = Pdo)
  }
} else {
  stop("Cannot find original data file, please check path.")
}

# --- 4. Custom labels and file mapping ---
# Note: list names are used for chart titles, file points to previously generated file names
metabolite_info <- list(
  Lactate = list(label = expression("Lactate (g/L)"), file = "Lactate_predictions.csv", col_name = "Lactate"),
  Acetate = list(label = expression("Acetate (g/L)"), file = "Acetate_predictions.csv", col_name = "Acetate"),
  Ethanol = list(label = expression("Ethanol (g/L)"), file = "Ethanol_predictions.csv", col_name = "Ethanol"),
  `1,3-PDO` = list(label = expression("1,3-PDO (g/L)"), file = "PDO_predictions.csv", col_name = "PDO")
)

# --- 5. Plot function (core logic) ---
plot_metabolite <- function(display_name, info) {
  
  # 5.1 Read prediction data
  pred_file <- file.path(result_dir, info$file)
  if(!file.exists(pred_file)) {
    warning(paste("Prediction file not found:", pred_file))
    return(NULL)
  }
  pred_df <- read_csv(pred_file, show_col_types = FALSE)
  
  # 5.2 Read real observed data
  # Use info$col_name to ensure correct column is found from real_data (e.g. "PDO")
  obs_df <- real_data %>% 
    select(time, value = !!sym(info$col_name)) %>% 
    filter(!is.na(value))
  
  # 5.3 Determine Training/Prediction boundary (take time midpoint)
  split_time <- median(obs_df$time)
  
  # 5.4 Dynamically calculate Y axis upper limit (ignore Frequentist extreme values)
  max_observed <- max(obs_df$value, na.rm = TRUE)
  max_bayes_ribbon <- max(pred_df$upper_bayes, na.rm = TRUE)
  
  # Take the larger of the two and add 30% space
  y_plot_limit <- max(max_observed, max_bayes_ribbon) * 1.3
  ymin <- 0 
  
  # 5.5 Plot
  p <- ggplot() +
    # (1) Background color block (Training area)
    annotate("rect", xmin = 0, xmax = split_time, ymin = -Inf, ymax = Inf,
             alpha = 0.2, fill = "grey85") +
    
    # (2) Frequentist (red)
    geom_ribbon(data = pred_df, aes(x = time, ymin = lower_nls, ymax = upper_nls),
                fill = "#FF0000", alpha = 0.2) + 
    geom_line(data = pred_df, aes(x = time, y = pred_nls), 
              color = "#FF0000", linewidth = 1, linetype = "dashed") +
    
    # (3) Bayesian (blue)
    geom_ribbon(data = pred_df, aes(x = time, ymin = lower_bayes, ymax = upper_bayes),
                fill = "#0000FF", alpha = 0.3) +
    geom_line(data = pred_df, aes(x = time, y = pred_bayes), 
              color = "#0000FF", linewidth = 1) +
    
    # (4) Real data points (black)
    geom_point(data = obs_df, aes(x = time, y = value), 
               size = 2.5, color = "black", shape = 19) +
    
    # (5) Text annotations
    annotate("text", x = split_time * 0.5, y = y_plot_limit * 0.95, 
             label = "Training", size = 4.5, fontface = "bold") +
    annotate("text", x = split_time + (75 - split_time) * 0.5, y = y_plot_limit * 0.95, 
             label = "Prediction", size = 4.5, fontface = "bold") +
    
    labs(x = "Time (h)", y = info$label) +
    
    # (6) Axis settings
    scale_x_continuous(
      name = "Time (h)",
      limits = c(0, 75),
      breaks = seq(0, 75, 15),
      expand = c(0, 0)
    ) +
    
    # (7) Force lock Y axis (clip lines that fly out)
    coord_cartesian(ylim = c(ymin, y_plot_limit), expand = FALSE) +
    
    theme_classic(base_size = 14) +
    theme(
      axis.line = element_line(linewidth = 0.8),
      axis.text = element_text(color = "black")
    )
  
  return(p)
}

# --- 6. Batch generate all plots ---
plots <- list()
for (met_name in names(metabolite_info)) {
  cat(paste("Generating plot for:", met_name, "...\n"))
  plots[[met_name]] <- plot_metabolite(met_name, metabolite_info[[met_name]])
}

# --- 7. Process metrics data (Bar Plot preparation) ---
# Read previously generated wide format data, convert to long format for ggplot
metrics_file <- file.path(result_dir, "model_metrics_comparison.csv")

if(file.exists(metrics_file)) {
  raw_metrics <- read_csv(metrics_file, show_col_types = FALSE)
  
  # Data transformation: Wide -> Long
  # We need to break RMSE_Frequentist, RMSE_Bayesian into Method and Metric
  metrics_df <- raw_metrics %>%
    pivot_longer(
      cols = -Metabolite, 
      names_to = c("Metric", "Method"), 
      names_sep = "_"
    ) %>%
    pivot_wider(names_from = Metric, values_from = value)
  
  # Process names (PDO -> 1,3-PDO) and set factor order
  metrics_df$Metabolite <- recode(metrics_df$Metabolite, "PDO" = "1,3-PDO")
  metrics_df$Metabolite <- factor(metrics_df$Metabolite, levels = c("Lactate", "Acetate", "Ethanol", "1,3-PDO"))
  
  # Generate Bar Plots (although not used in main plot, keep logic)
  bar_rmse <- ggplot(metrics_df, aes(x = Metabolite, y = RMSE, fill = Method)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("Frequentist" = "#FF0000", "Bayesian" = "#0000FF")) +
    theme_minimal(base_size = 12) + labs(y = "RMSE") + theme(legend.position = "none")
  
  bar_nsme <- ggplot(metrics_df, aes(x = Metabolite, y = NSME, fill = Method)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("Frequentist" = "#FF0000", "Bayesian" = "#0000FF")) +
    theme_minimal(base_size = 12) + labs(y = "NSME") + theme(legend.position = "none")
}

# --- 8. Generate legend ---
legend_plot <- ggplot() +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1), fill = NA) +
  # Frequentist
  annotate("rect", xmin = 0.15, xmax = 0.20, ymin = 0.4, ymax = 0.6, fill = "#FF0000", alpha=0.3) +
  annotate("segment", x = 0.15, xend = 0.20, y = 0.5, yend = 0.5, color = "#FF0000", linetype="dashed", linewidth=1) +
  annotate("text", x = 0.22, y = 0.5, label = "Frequentist (NLS)", hjust = 0, size = 5) +
  # Bayesian
  annotate("rect", xmin = 0.55, xmax = 0.60, ymin = 0.4, ymax = 0.6, fill = "#0000FF", alpha=0.3) +
  annotate("segment", x = 0.55, xend = 0.60, y = 0.5, yend = 0.5, color = "#0000FF", linewidth=1) +
  annotate("text", x = 0.62, y = 0.5, label = "Bayesian Inference", hjust = 0, size = 5) +
  theme_void()

# --- 9. Final assembly and save ---
# Check if all plots are generated successfully
valid_plots <- plots[!sapply(plots, is.null)]

if (length(valid_plots) == 4) {
  main_plot <- plot_grid(
    plots$Lactate, plots$Acetate,
    plots$Ethanol, plots$`1,3-PDO`,
    labels = c("A", "B", "C", "D"),
    ncol = 2, label_size = 16, align = "hv"
  )
  
  final_plot <- plot_grid(
    main_plot,
    legend_plot,
    ncol = 1, rel_heights = c(1, 0.08) # 调整底部图例高度
  )
  
  # 保存
  save_path <- file.path(result_dir, "Figure5_Prediction_Comparison.png")
  ggsave(save_path, final_plot, width = 12, height = 12, type = "cairo",dpi = 1000)
  
  print(paste("可视化完成！文件已保存至:", save_path))
  print(final_plot) # 在 RStudio 中显示
} else {
  warning("未能生成所有代谢物的图，请检查是否有部分代谢物预测失败。")
}