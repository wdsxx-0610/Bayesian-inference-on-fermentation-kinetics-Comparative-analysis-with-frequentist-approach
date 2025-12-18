graphics.off()
rm(list = ls())

library(tidyverse)
library(scales)
library(patchwork)
library(cowplot)

# =============================
# Basic parameters
# =============================
fileName <- "./results/"
cofermentation_names <- c("Lactate", "Acetate", "Ethanol", "PDO")

# =============================
# Plot function (with scale control & automatic expansion)
# =============================
generate_plot <- function(metabolite, y_limits = NULL, y_breaks = NULL) {
  
  # Read data
  df_f <- read.csv(file.path(fileName, metabolite, paste0(metabolite, "_frequentist.csv")))
  df_b <- read.csv(file.path(fileName, metabolite, paste0(metabolite, "_1211bayesian.csv")))
  
  pred_col_f <- paste0(metabolite, ".1")
  
  obsData <- tibble(
    Time = df_f$time[1:28],
    Value = df_f[[metabolite]][1:28]
  )
  
  data_f <- tibble(
    Time = df_f$time.1,
    Value = df_f[[pred_col_f]],
    Lower = df_f$lwr,
    Upper = df_f$upr,
    Method = "Frequentist"
  )
  
  data_b <- tibble(
    Time = df_b$time.1,
    Value = df_b$estimate,
    Lower = df_b$conf.low,
    Upper = df_b$conf.high,
    Method = "Bayesian"
  )
  
  combined_data <- bind_rows(data_f, data_b)
  
  # If no y_limits given, calculate automatically
  if (is.null(y_limits)) {
    ymin_auto <- min(combined_data$Lower, na.rm = TRUE)
    ymax_auto <- max(combined_data$Upper, na.rm = TRUE) * 1.05
    y_limits <- c(ymin_auto, ymax_auto)
  }
  
  # If y_limits manually given but data max exceeds upper limit, expand automatically
  max_data <- max(combined_data$Upper, na.rm = TRUE)
  if (y_limits[2] < max_data) {
    y_limits[2] <- max_data * 1.05
  }
  
  # Automatically calculate scale
  if (is.null(y_breaks)) {
    y_breaks <- pretty(y_limits, n = 5)
  }
  
  # Set y-axis label
  label_text <- ifelse(metabolite == "PDO",
                       '"1,3-PDO"~(g/L)',
                       paste0('"', metabolite, '"~(g/L)'))
  # --- x-axis (fixed max 75) ---
  x_limits <- c(0, 78)
  x_breaks <- seq(0, 75, 15)
  
  # Plot
  p <- ggplot() +
    geom_ribbon(data = combined_data,
                aes(x = Time, ymin = Lower, ymax = Upper, fill = Method),
                alpha = 0.3) +
    geom_line(data = combined_data,
              aes(x = Time, y = Value, color = Method), size = 0.5) +
    geom_point(data = obsData,
               aes(x = Time, y = Value), shape = 21, fill = "black", size = 1) +
    scale_fill_manual(values = c("Frequentist" = "#FF0000", "Bayesian" = "#0000FF")) +
    scale_color_manual(values = c("Frequentist" = "#FF0000", "Bayesian" = "#0000FF")) +
    scale_x_continuous(name = "Time (h)",
                       limits = x_limits,
                       breaks = x_breaks,
                       expand = c(0, 0)) +
    scale_y_continuous(
      name = parse(text = label_text),
      limits = y_limits,
      breaks = y_breaks,
      expand = c(0, 0),
      labels = scales::label_number(accuracy = 0.1)
    ) +
    facet_wrap(~Method, ncol = 2) +
    theme_classic(base_size = 13) +
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      legend.position = "none"
    )
  
  return(p)
}

# =============================
# Manual scale settings
# =============================
y_settings <- list(
  Lactate = list(limits = c(0, 16.0), breaks = seq(0, 16.0, 4.0)),
  Acetate = list(limits = c(0, 60.0), breaks = seq(0, 60.0, 15.0)),
  Ethanol = list(limits = c(0, 8.0), breaks = seq(0, 8.0, 2.0)),
  PDO     = list(limits = c(0, 90.0), breaks = seq(0, 90.0, 15.0))
)

# Generate main plot list
plot_list <- lapply(cofermentation_names, function(met) {
  generate_plot(met,
                y_limits = y_settings[[met]]$limits,
                y_breaks = y_settings[[met]]$breaks)
})

# =============================
# Read model evaluation metrics
# =============================
metrics_df <- data.frame(Metabolite = character(),
                         Method = character(),
                         RMSE = numeric(),
                         NSME = numeric(),
                         stringsAsFactors = FALSE)

for (metabolite in cofermentation_names) {
  bayes_summary <- read.csv(paste0(fileName, metabolite, "/", metabolite, "_bayesian_summary.csv"))
  freq_summary  <- read.csv(paste0(fileName, metabolite, "/", metabolite, "_frequentist_summary.csv"))
  
  rmse_b <- as.numeric(bayes_summary[8, 2])
  nsme_b <- as.numeric(bayes_summary[10, 2])
  rmse_f <- as.numeric(freq_summary[6, 2])
  nsme_f <- as.numeric(freq_summary[8, 2])
  
  metrics_df <- metrics_df %>%
    add_row(Metabolite = metabolite, Method = "Bayesian", RMSE = rmse_b, NSME = nsme_b) %>%
    add_row(Metabolite = metabolite, Method = "Frequentist", RMSE = rmse_f, NSME = nsme_f)
}

metrics_long <- metrics_df %>%
  pivot_longer(cols = c("RMSE", "NSME"), names_to = "Metric", values_to = "Value") %>%
  mutate(
    Method = factor(Method, levels = c("Frequentist", "Bayesian")),
    Metric = factor(Metric, levels = c("RMSE", "NSME")),
    Metabolite = recode(Metabolite, PDO = "1,3-PDO"),
    Metabolite = factor(Metabolite, levels = c("Lactate", "Acetate", "Ethanol", "1,3-PDO"))
  )

# =============================
# Bar charts
# =============================
bar_rmse <- ggplot(filter(metrics_long, Metric == "RMSE"),
                   aes(x = Metabolite, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.6) +
  scale_fill_manual(values = c("Frequentist" = "#FF0000", "Bayesian" = "#0000FF")) +
  scale_y_continuous(limits = c(0, 1.75), breaks = seq(0, 1.75, 0.25)) +  
  theme_minimal(base_size = 13) +
  labs(y = "RMSE", x = NULL, fill = "Method") +
  theme(legend.position = "none")

bar_nsme <- ggplot(filter(metrics_long, Metric == "NSME"),
                   aes(x = Metabolite, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.6) +
  scale_fill_manual(values = c("Frequentist" = "#FF0000", "Bayesian" = "#0000FF")) +
  theme_minimal(base_size = 13) +
  labs(y = "NSME", x = NULL, fill = "Method") +
  theme(legend.position = "none")

legend_plot <- ggplot() +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1), fill = NA) +
  annotate("rect", xmin = 0.05, xmax = 0.12, ymin = 0.4, ymax = 0.6, fill = "#FF0000") +
  annotate("text", x = 0.14, y = 0.5, label = "Frequentist", hjust = 0, size = 3.5) +
  annotate("rect", xmin = 0.35, xmax = 0.42, ymin = 0.4, ymax = 0.6, fill = "#0000FF") +
  annotate("text", x = 0.44, y = 0.5, label = "Bayesian", hjust = 0, size = 3.5) +
  theme_void()



library(ggplot2)
library(cowplot)

# --- 1. Define single-row layout parameters ---
y_pos <- 0.5       # All elements on the same horizontal line
line_len <- 0.05   # Width of lines/blocks
gap_text <- 0.01   # Gap between icon and text

# Colors
col_bayes <- "#0000FF"
col_freq  <- "#FF0000"

legend_G_plot <- ggplot() +
  # =======================================================
# Group 1: Bayesian blue line (Bayesian Mean) - Position X ~ 0.05
# =======================================================
annotate("segment", x = 0.05, xend = 0.05 + line_len, y = y_pos, yend = y_pos, 
         color = col_bayes, linewidth = 1.2) +
  annotate("text", x = 0.05 + line_len + gap_text, y = y_pos, 
           label = "Bayesian estimate mean", 
           hjust = 0, vjust = 0.5, size = 4, family = "Arial", color = "#333333") +
  
  # =======================================================
# Group 2: Bayesian blue box (Bayesian CI) - Position X ~ 0.28
# =======================================================
annotate("rect", xmin = 0.28, xmax = 0.28 + line_len, 
         ymin = y_pos - 0.15, ymax = y_pos + 0.15, 
         fill = col_bayes, alpha = 0.3) +
  annotate("text", x = 0.28 + line_len + gap_text, y = y_pos, 
           label = "95% high density interval", 
           hjust = 0, vjust = 0.5, size = 4, family = "Arial", color = "#333333") +
  
  # =======================================================
# Group 3: Frequentist red line (Freq Mean) - Position X ~ 0.53
# =======================================================
annotate("segment", x = 0.53, xend = 0.53 + line_len, y = y_pos, yend = y_pos, 
         color = col_freq, linewidth = 1.2) +
  annotate("text", x = 0.53 + line_len + gap_text, y = y_pos, 
           label = "Frequentist estimate", 
           hjust = 0, vjust = 0.5, size = 4, family = "Arial", color = "#333333") +
  
  # =======================================================
# Group 4: Frequentist red box (Freq CI) - Position X ~ 0.75
# =======================================================
annotate("rect", xmin = 0.75, xmax = 0.75 + line_len, 
         ymin = y_pos - 0.15, ymax = y_pos + 0.15, 
         fill = col_freq, alpha = 0.3) +
  annotate("text", x = 0.75 + line_len + gap_text, y = y_pos, 
           label = "95% confidence interval", 
           hjust = 0, vjust = 0.5, size = 4, family = "Arial", color = "#333333") +
  
  # --- Set canvas ---
  scale_x_continuous(limits = c(0, 1), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  theme_void()
# 1. Generate grid for upper 4 plots (A, B, C, D)
main_plots <- plot_grid(
  plot_list[[1]], plot_list[[2]],
  plot_list[[3]], plot_list[[4]],
  labels = LETTERS[1:4], 
  label_size = 14,
  ncol = 2, 
  align = 'hv' # Ensure alignment
)

# 2. Vertically combine main plots and legend
final_plot <- plot_grid(
  main_plots,     # Upper part
  legend_G_plot,  # Lower part
  ncol = 1,
  # Adjust height ratio:
  # Upper plots take most space (e.g. 20), legend takes small part (e.g. 1)
  rel_heights = c(20, 1) 
)

print(final_plot)

# 3. 保存
# 注意：因为去掉了柱状图，高度(height)可能需要稍微减小一点，避免留白太多
ggsave(filename = paste0(fileName, "Figure2_final_4panel.png"),
       plot = final_plot, type = "cairo",
       width = 12, height = 10, dpi = 1000)
# 保
