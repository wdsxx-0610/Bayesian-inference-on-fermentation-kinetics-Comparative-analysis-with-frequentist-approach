# ==============================================================================
# Visualization: Frequentist vs Bayesian Predictions
# ==============================================================================

# --- 1. 加载包 ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, cowplot, scales, readr)

# --- 2. 设置路径 (匹配之前的运行结果) ---
# 预测结果所在的文件夹
result_dir <- "/Users/wdsxx0610/Documents/R_directory/A_formation/cofermentation/Half_prediction_Results"
# 原始数据文件路径
raw_data_path <- "/Users/wdsxx0610/Documents/R_directory/A_formation/cofermentation/CofermentationData.csv"

# --- 3. 读取原始数据 ---
if(file.exists(raw_data_path)){
  real_data <- read.csv(raw_data_path) %>%
    rename(time = time) %>% # 确保时间列叫 time
    rename_with(~ str_to_title(.), .cols = -time) # 确保列名大写 (Lactate, etc.)
  
  # 修复 PDO 列名 (如果变成了 Pdo -> PDO)
  if ("Pdo" %in% colnames(real_data)) {
    real_data <- real_data %>% rename(PDO = Pdo)
  }
} else {
  stop("无法找到原始数据文件，请检查路径。")
}

# --- 4. 自定义标签与文件映射 ---
# 注意：列表的名字用于图表标题，file 指向之前生成的文件名
metabolite_info <- list(
  Lactate = list(label = expression("Lactate (g/L)"), file = "Lactate_predictions.csv", col_name = "Lactate"),
  Acetate = list(label = expression("Acetate (g/L)"), file = "Acetate_predictions.csv", col_name = "Acetate"),
  Ethanol = list(label = expression("Ethanol (g/L)"), file = "Ethanol_predictions.csv", col_name = "Ethanol"),
  `1,3-PDO` = list(label = expression("1,3-PDO (g/L)"), file = "PDO_predictions.csv", col_name = "PDO")
)

# --- 5. 绘图函数 (核心逻辑) ---
plot_metabolite <- function(display_name, info) {
  
  # 5.1 读取预测数据
  pred_file <- file.path(result_dir, info$file)
  if(!file.exists(pred_file)) {
    warning(paste("Prediction file not found:", pred_file))
    return(NULL)
  }
  pred_df <- read_csv(pred_file, show_col_types = FALSE)
  
  # 5.2 读取真实观测数据
  # 使用 info$col_name 来确保能从 real_data 里找到正确的列 (比如 "PDO")
  obs_df <- real_data %>% 
    select(time, value = !!sym(info$col_name)) %>% 
    filter(!is.na(value))
  
  # 5.3 确定 Training/Prediction 分界线 (取时间中点)
  split_time <- median(obs_df$time)
  
  # 5.4 动态计算 Y 轴上限 (忽略 Frequentist 的极端值)
  max_observed <- max(obs_df$value, na.rm = TRUE)
  max_bayes_ribbon <- max(pred_df$upper_bayes, na.rm = TRUE)
  
  # 取两者较大值，并增加 30% 空间
  y_plot_limit <- max(max_observed, max_bayes_ribbon) * 1.3
  ymin <- 0 
  
  # 5.5 绘图
  p <- ggplot() +
    # (1) 背景色块 (Training区域)
    annotate("rect", xmin = 0, xmax = split_time, ymin = -Inf, ymax = Inf,
             alpha = 0.2, fill = "grey85") +
    
    # (2) Frequentist (红色)
    geom_ribbon(data = pred_df, aes(x = time, ymin = lower_nls, ymax = upper_nls),
                fill = "#FF0000", alpha = 0.2) + 
    geom_line(data = pred_df, aes(x = time, y = pred_nls), 
              color = "#FF0000", linewidth = 1, linetype = "dashed") +
    
    # (3) Bayesian (蓝色)
    geom_ribbon(data = pred_df, aes(x = time, ymin = lower_bayes, ymax = upper_bayes),
                fill = "#0000FF", alpha = 0.3) +
    geom_line(data = pred_df, aes(x = time, y = pred_bayes), 
              color = "#0000FF", linewidth = 1) +
    
    # (4) 真实数据点 (黑色)
    geom_point(data = obs_df, aes(x = time, y = value), 
               size = 2.5, color = "black", shape = 19) +
    
    # (5) 文字标注
    annotate("text", x = split_time * 0.5, y = y_plot_limit * 0.95, 
             label = "Training", size = 4.5, fontface = "bold") +
    annotate("text", x = split_time + (75 - split_time) * 0.5, y = y_plot_limit * 0.95, 
             label = "Prediction", size = 4.5, fontface = "bold") +
    
    labs(x = "Time (h)", y = info$label) +
    
    # (6) 坐标轴设置
    scale_x_continuous(
      name = "Time (h)",
      limits = c(0, 75),
      breaks = seq(0, 75, 15),
      expand = c(0, 0)
    ) +
    
    # (7) 强制锁定 Y 轴 (剪裁掉飞出去的线)
    coord_cartesian(ylim = c(ymin, y_plot_limit), expand = FALSE) +
    
    theme_classic(base_size = 14) +
    theme(
      axis.line = element_line(linewidth = 0.8),
      axis.text = element_text(color = "black")
    )
  
  return(p)
}

# --- 6. 批量生成所有图片 ---
plots <- list()
for (met_name in names(metabolite_info)) {
  cat(paste("Generating plot for:", met_name, "...\n"))
  plots[[met_name]] <- plot_metabolite(met_name, metabolite_info[[met_name]])
}

# --- 7. 处理指标数据 (Bar Plot 准备) ---
# 读取之前生成的宽格式数据，转换为长格式以适配 ggplot
metrics_file <- file.path(result_dir, "model_metrics_comparison.csv")

if(file.exists(metrics_file)) {
  raw_metrics <- read_csv(metrics_file, show_col_types = FALSE)
  
  # 数据转换: Wide -> Long
  # 我们需要把 RMSE_Frequentist, RMSE_Bayesian 拆解成 Method 和 Metric
  metrics_df <- raw_metrics %>%
    pivot_longer(
      cols = -Metabolite, 
      names_to = c("Metric", "Method"), 
      names_sep = "_"
    ) %>%
    pivot_wider(names_from = Metric, values_from = value)
  
  # 处理名称 (PDO -> 1,3-PDO) 并设置因子顺序
  metrics_df$Metabolite <- recode(metrics_df$Metabolite, "PDO" = "1,3-PDO")
  metrics_df$Metabolite <- factor(metrics_df$Metabolite, levels = c("Lactate", "Acetate", "Ethanol", "1,3-PDO"))
  
  # 生成 Bar Plots (虽然主图没用到，但保留逻辑)
  bar_rmse <- ggplot(metrics_df, aes(x = Metabolite, y = RMSE, fill = Method)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("Frequentist" = "#FF0000", "Bayesian" = "#0000FF")) +
    theme_minimal(base_size = 12) + labs(y = "RMSE") + theme(legend.position = "none")
  
  bar_nsme <- ggplot(metrics_df, aes(x = Metabolite, y = NSME, fill = Method)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("Frequentist" = "#FF0000", "Bayesian" = "#0000FF")) +
    theme_minimal(base_size = 12) + labs(y = "NSME") + theme(legend.position = "none")
}

# --- 8. 生成图例 (Legend) ---
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

# --- 9. 最终拼图与保存 ---
# 检查是否所有图都生成成功
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