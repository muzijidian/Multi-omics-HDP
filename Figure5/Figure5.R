rm(list = ls())
setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code')
library(dplyr)
library(stringr)
library(reshape2)
library(openxlsx)
library(ggplot2)
library(ggh4x)
library(patchwork)
source('help_func.R')


omic_level = c("Baseline", 'Metabolome', 'Fungi', 'Bacteria', 'Fun+Bac',
               'Met+Fun',  'Met+Bac', 
               'Multi-Omics')

type_col = c(
  "Baseline" = "#bebebe",
  'Metabolome' = "#4A7B9D", 
  'Fungi' = "#6E7582",
  'Bacteria' = "#839B97",
  'Fun+Bac' = "#B6A39E",
  'Met+Fun' = "#9B6B6C",
  'Met+Bac' = "#E69F00",
  'Multi-Omics' = "#B24745"
)

col_outcome = c('PE'='#e3a9a4', 'GH'='#f0d575', 'PE/GH'='#edad8a')

# 定义一个【万能清洗函数】：不管名字有没有变异，强制改回标准格式
fix_names <- function(x) {
  x <- gsub("\\.", "+", x)       # 把所有的点 . 变回 +
  x <- gsub("Multi\\+Omics", "Multi-Omics", x) # 修正 Multi-Omics
  return(x)
}

# 定义一个【安全映射函数】：不管Excel里是旧名字还是新名字，都统一
safe_map_outcome <- function(x) {
  case_when(
    x %in% c('GH_vs_NP', 'GH') ~ 'GH',
    x %in% c('PE_vs_NP', 'PE') ~ 'PE',
    x %in% c('Combined_vs_NP', 'PE/GH') ~ 'PE/GH',
    TRUE ~ x
  )
}

# 2. 绘制 ROC 曲线 (左图) =====================================================
roc_file <- 'Anti_results/prediction_roc_curve_data_upsampling_best.xlsx'
roc_data <- read.xlsx(roc_file)

# 【修正步骤 1】安全转换 Outcome
roc_data$outcome <- safe_map_outcome(roc_data$outcome)

# 指定显示顺序
roc_data$outcome <- factor(roc_data$outcome, levels = c('PE', 'GH', 'PE/GH'))

# 【修正步骤 2】清洗 ft_type 并设置因子
roc_data$ft_type <- fix_names(roc_data$ft_type)
roc_data$ft_type <- factor(roc_data$ft_type, levels = omic_level)

# 调试检查：确保这里不是0
print(table(roc_data$outcome))

# 绘图
p1 <- ggplot(roc_data, aes(x = fpr, y = mean_sens, color = ft_type)) +
  facet_wrap2(~outcome, 
              strip = strip_themed(background_x = elem_list_rect(fill = col_outcome))) +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey30") +
  scale_color_manual(values = type_col) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "1 - Specificity", y = "Sensitivity", color = "Models") +
  # coord_equal() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18, color = "black"), 
    axis.text = element_text(size = 16, color = "grey10"),
    panel.grid = element_line(linetype = 3, color = 'grey95'),
    # legend.position = 'bottom',
    legend.position = "none",
    legend.text = element_text(size = 8),
    strip.text = element_text(size = 18, color = 'white', face = "bold"), 
    panel.spacing = unit(5, "mm"),
    plot.margin = ggplot2::margin(5, 5, 5, 5)
  )

# 3. 绘制 AUC 热图 (右图) =====================================================
auc_file <- 'Anti_results/prediction_summary_auc_upsampling_best.xlsx'
auc_data <- read.xlsx(auc_file)

# 【修正步骤 3】对 Heatmap 数据做同样的安全处理
auc_data$outcome <- safe_map_outcome(auc_data$outcome)

auc_data$outcome <- factor(auc_data$outcome, levels = c('PE', 'GH', 'PE/GH'))

auc_data$ft_type <- fix_names(auc_data$ft_type)
auc_data$ft_type <- factor(auc_data$ft_type, levels = omic_level)

auc_data <- auc_data %>% 
  mutate(auc_label = sprintf('%.3f\n(%.3f)', roc_auc, auc_sd))
y_axis_colors <- unname(type_col[rev(omic_level)])
# 绘图
p2 <- ggplot(auc_data, aes(x = outcome, y = ft_type)) +
  geom_tile(aes(fill = roc_auc), color = 'white') +
  geom_text(aes(label = auc_label, color = ifelse(roc_auc > 0.8, 'white', 'black')), 
            show.legend = FALSE, size = 3.8) +
  scale_color_manual(values = c('black', 'white')) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(omic_level), expand = c(0, 0)) + 
  scale_fill_gradientn(
    colors = c("#F2F7FA", "#D4E4ED", "#AFC9D9", "#7DA8C3", "#4A7B9D", "#28506B", "#0E2F40"),
    limits = c(0.5, 1.0), 
    breaks = seq(0.5, 1, 0.1),
    oob = scales::squish,
    name = "AUC (SD)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14, color = 'black', angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 11, color = y_axis_colors, face = "bold"),
    # axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.6, "cm"),
    legend.key.width = unit(0.6, "cm"),
    legend.position = "right",
    plot.margin = ggplot2::margin(5, 5, 5, 0)
  )

# 4. 拼图并保存 ===============================================================
combined_plot <- (p1 | p2) + plot_layout(widths = c(3, 0.7))

print(combined_plot)


ggsave('Anti_graph/Fig_New_Analysis_ROC_AUC_upsampling_best_1.pdf', combined_plot, width = 15, height = 4.2)













rm(list = ls())
setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code')
library(ggplot2)
library(dplyr)
library(readxl)
library(openxlsx)
library(stringr)
library(ggh4x) # 核心包：用于 facet_manual 和 strip_themed


# 重新读取之前的特征列表，以便给 Feature 打标签
info_met <- read_excel('Anti_results/Met_sig.xlsx')
info_its <- read_excel('Anti_results/ITS_sig.xlsx', sheet = 1)
info_bac <- read_excel('Anti_results/Taxonomy_s_sig.xlsx', sheet = 1)

# 获取 ID 列表用于分类 source (保持原有逻辑)
f_met <- info_met[['meta_name']]
f_its <- info_its[['meta_name']]
f_bac <- info_bac[['meta_name']]
cov_names <- c("age", "parity", "edu", "smk", "drk", "wk", "BMI_prep") 

# 读取你的重要性结果
imp_tab_raw <- read.xlsx('Anti_results/prediction_importance_upsampling_best.xlsx')

imp_tab_raw <- imp_tab_raw %>%
  mutate(Feature = case_when(
    Feature == "parity1" ~ "parity",
    Feature == "drk1"    ~ "drk",
    Feature == "edu1"    ~ "edu",
    Feature == "smk1"    ~ "smk",
    TRUE                 ~ Feature # 其他变量保持原样
  ))

# ==============================================================================
# 2. 数据清洗与整理
# ==============================================================================

# 2.1 筛选 Multi-Omics 模型，并添加 Source 列
plot_data <- imp_tab_raw %>%
  filter(ft_type == 'Multi-Omics') %>%  # 只画多组学整合模型
  mutate(
    source = case_when(
      Feature %in% f_met ~ "Metabolites",
      Feature %in% f_its ~ "Fungi",
      Feature %in% f_bac ~ "Bacteria",
      Feature %in% cov_names ~ "Covariates",
      TRUE ~ "Other" # 防止有漏网之鱼
    )
  )


plot_data <- plot_data %>%
  # 1. 匹配代谢物名称
  left_join(info_met %>% select(meta_name, sugg_cmpd_name), by = c("Feature" = "meta_name")) %>%
  # 2. 匹配真菌名称
  left_join(info_its %>% select(meta_name, name_its = name), by = c("Feature" = "meta_name")) %>%
  # 3. 匹配细菌名称
  left_join(info_bac %>% select(meta_name, name_bac = name), by = c("Feature" = "meta_name")) %>%
  mutate(
    # 优先使用匹配到的真实名称，如果没有匹配到(比如Baseline)则用Feature ID
    real_name_raw = case_when(
      source == "Metabolites" ~ sugg_cmpd_name,
      source == "Fungi"      ~ name_its,
      source == "Bacteria"   ~ name_bac,
      TRUE                   ~ Feature 
    ),
    # 如果对应的列是 NA (有的可能没注释)，回退到原始 ID
    real_name_raw = coalesce(real_name_raw, Feature)
  )



specific_names_map <- c(
  '2,4-diacetamino-2,4,6-triphenoxy-D-mannopyranose' = 'Man2NAc4NAc',
  'Glycerophospho-N-Arachidonoyl Ethanolamine' = 'GP-NAE',
  'Phosphatidylethanolamine lyso alkenyl 16:0' = 'LPE(P-16:0)',
  'Glycine deoxycholic acid'= 'GDCA',
  'Glycochenodeoxycholic acid'= 'GCDCA',
  'wk' = "Gestational age",
  'age' = "Maternal age",
  'BMI_prep' = "Pre-pregnancy BMI",
  'smk' = "Smoking status",
  'drk' = "Drinking status",
  'parity' = "Parity",
  'edu' = "Education level"
  
)

plot_data <- plot_data %>%
  mutate(
    name_display = coalesce(specific_names_map[real_name_raw], real_name_raw),

    name_display = str_replace_all(name_display, c(" Acid" = " acid", " ACID" = " acid")),

    name_display = name_display %>%
      str_remove_all("^s__|^g__|^p__|^f__") %>% 
      str_replace_all("_", " ") %>%
      str_trunc(35) 
  )




# 2.3 计算归一化重要性 (Normalized Importance)
# 逻辑：在每个 Outcome (GH, PE, Combined) 内部，将 importance 归一化为比例
plot_data <- plot_data %>%
  group_by(outcome) %>%
  mutate(imp_value = Overall / sum(Overall)) %>% # 计算相对比例
  ungroup()

# 2.4 筛选 Top 特征并排序
# 逻辑：计算每个特征在所有 Outcome 中的重要性总和，取 Top 30 (或40，可视情况调整)
top_features <- plot_data %>%
  group_by(name_display, source) %>%
  summarise(total_imp = sum(imp_value), .groups = 'drop') %>%
  arrange(desc(total_imp)) 
# %>%
#   slice_head(n = 40) # *** 这里修改显示的特征数量 ***


# top_features_stats <- plot_data %>%
#   group_by(Feature, name_display, source) %>%
#   summarise(total_imp = sum(imp_value), .groups = 'drop') %>%
#   arrange(desc(total_imp))
# 
# met_top10_ids <- top_features_stats %>%
#   filter(source == "Metabolome") %>%
#   arrange(desc(total_imp)) %>%
#   slice_head(n = 10) %>%
#   pull(Feature) # 获取原始 ID (meta_name)
# 
# # 从原始注释表中提取这10个代谢物的完整信息
# met_sig_top10_df <- info_met %>%
#   filter(meta_name %in% met_top10_ids) %>%
#   # 按照重要性顺序排列 (可选)
#   mutate(meta_name = factor(meta_name, levels = met_top10_ids)) %>%
#   arrange(meta_name)
# 
# # 保存文件
# write.xlsx(met_sig_top10_df, 'Anti_results/Met_sig_top10.xlsx')


plot_data_final <- plot_data %>%
  inner_join(top_features %>% select(name_display), by = "name_display")

source_levels <- c('Metabolites', 'Fungi', 'Bacteria', 'Covariates')

plot_data_final <- plot_data_final %>%
  mutate(outcome = case_when(
    outcome %in% c('PE_vs_NP', 'PE') ~ 'PE',
    outcome %in% c('GH_vs_NP', 'GH') ~ 'GH',
    outcome %in% c('Combined_vs_NP', 'PE/GH') ~ 'PE/GH',
    TRUE ~ as.character(outcome)
  ))

# 2. 设置特征 (Y轴) 的顺序 (保持不变)
ft_order <- top_features %>% pull(name_display)
plot_data_final$name_display <- factor(plot_data_final$name_display, levels = ft_order)

plot_data_final$source <- factor(plot_data_final$source, levels = source_levels)

plot_data_final$outcome <- factor(plot_data_final$outcome, 
                                  levels = c( "PE", "GH","PE/GH"))




# ==============================================================================
# 3. 动态计算分面宽度
# ==============================================================================

library(tidyr)
n_counts <- plot_data_final %>% 
  select(name_display, source) %>% 
  distinct() %>% 
  count(source) %>%
  mutate(source = factor(source, levels = source_levels)) %>% 
  arrange(source) %>%
  complete(source, fill = list(n = 1)) 

facet_widths <- n_counts$n + 2 
design_str <- paste(LETTERS[1:length(source_levels)], collapse = "")

# ==============================================================================
# 4. 绘图
# ==============================================================================

col_outcome <- c('PE'='#e3a9a4', 'GH'='#f6d364', 'PE/GH'='#edad8a')
col_source_map <- c(
  'Metabolites' = "#4A7B9D",
  'Fungi'      = "#f0d9a2", 
  'Bacteria'   = "#ef8e86",
  'Covariates'   = "grey70"
)
col_source_vec <- col_source_map[levels(plot_data_final$source)]

p <- ggplot(plot_data_final, aes(x = name_display, y = imp_value, fill = outcome)) +
  facet_manual(
    ~source, 
    scales = "free_x", 
    design = design_str, 
    widths = facet_widths, 
    strip = strip_themed(
      background_x = elem_list_rect(fill = col_source_vec)
    )
  ) +
  geom_col(position = position_stack(), width = 0.79) +
  scale_fill_manual(values = col_outcome, name = "") +
  scale_y_continuous(expand = c(0, 0)) + 
  labs(y = "Normalized importance") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    # 调整文字角度和大小以适应较长的真实名称
    axis.text.x = element_text(size = 12, color = "black", angle = -40, vjust = 1, hjust = 0), 
    
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 14, color = "grey10", vjust = 0.5, hjust = 0.5),
    
    axis.ticks.length = unit(1.5, 'mm'), 
    panel.grid = element_blank(), 
    
    strip.text = element_text(size = 14, color = 'white', face = "bold"), 
    strip.placement = 'inside',
    panel.spacing = unit(2, "mm"), 
    
    legend.position = 'right',
    legend.background = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 12),
    
    plot.margin = ggplot2::margin(5, 5, 5, 5, 'mm')
  )

print(p)

# 保存 (增加了高度，因为名字旋转后占空间更大)
ggsave('Anti_graph/Fig_Importance_Stacked_Names.pdf', p, width = 16, height = 5.5)








rm(list = ls())
setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code') # 请确保路径正确

library(dplyr)
library(ggplot2)
library(readxl)
library(patchwork)
library(stringr)

# ==============================================================================
# 1. 全局设置：颜色与模型顺序
# ==============================================================================

# 指定模型顺序
target_models <- c("Baseline", "Met+Fun", "Met+Bac", "Multi-Omics")

# 模型线条颜色
type_col <- c(
  "Baseline" = "#bebebe",
  "Met+Fun" = "#9B6B6C",
  "Met+Bac" = "#E69F00",
  "Multi-Omics" = "#B24745"
)

# 定义通用的分面文字样式
strip_text_style <- element_text(size = 16, color = 'white', face = "bold")

# 定义通用主题 (去除背景网格，保留边框)
common_theme <- theme_bw() +
  theme(
    axis.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    strip.text = strip_text_style 
  )



# ==============================================================================
# 2. 处理 THSBC 数据 (左图 - 保留 Y 轴)
# ==============================================================================

data_thsbc <- read_excel('Anti_results/prediction_roc_curve_data_THSBC.xlsx')

plot_data_thsbc <- data_thsbc %>%
  filter(outcome %in% c('PE_vs_NP', 'PE')) %>%
  filter(ft_type %in% target_models) %>%
  mutate(ft_type = factor(ft_type, levels = target_models)) %>%
  rename(sensitivity = mean_sens) %>%
  mutate(Cohort = "THSBC (PE)")

# 绘图 P1 (保留完整 Y 轴)
p1 <- ggplot(plot_data_thsbc, aes(x = fpr, y = sensitivity, color = ft_type)) +
  facet_wrap(~Cohort) + 
  geom_line(linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = type_col) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  common_theme +
  theme(strip.background = element_rect(fill = "#d66a43", color = NA))

# ==============================================================================
# 3. 处理 WeBirth 数据 (中图 - 去除 Y 轴)
# ==============================================================================

data_webirth <- read.csv('Anti_results/Final_ROC_Curves_Data_WeBirth.csv')

if("sens" %in% colnames(data_webirth)) {
  data_webirth <- data_webirth %>% rename(mean_sens = sens)
}

plot_data_webirth <- data_webirth %>%
  filter(outcome %in% c('PE_vs_NP', 'PE')) %>%
  filter(ft_type %in% target_models) %>%
  mutate(ft_type = factor(ft_type, levels = target_models)) %>%
  rename(sensitivity = mean_sens) %>%
  mutate(Cohort = "WeBirth (PE)")

# 绘图 P2
p2 <- ggplot(plot_data_webirth, aes(x = fpr, y = sensitivity, color = ft_type)) +
  facet_wrap(~Cohort) +
  geom_line(linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = type_col) +
  # [修改 1] y = NULL 去掉标题
  labs(x = "1 - Specificity", y = NULL) + 
  common_theme +
  theme(
    strip.background = element_rect(fill = "#a6bfbd", color = NA),
    # [修改 2] 去掉 Y 轴刻度文字和刻度线
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )

# ==============================================================================
# 4. 处理 HBC 数据 (右图 - 去除 Y 轴)
# ==============================================================================

data_HBC <- read.csv('Anti_results/Final_ROC_Curves_Data_HBC.csv')

if("sens" %in% colnames(data_HBC)) {
  data_HBC <- data_HBC %>% rename(mean_sens = sens)
}

data_HBC <- data_HBC %>%
  filter(outcome %in% c('PE_vs_NP', 'PE')) %>%
  filter(ft_type %in% target_models) %>%
  mutate(ft_type = factor(ft_type, levels = target_models)) %>%
  rename(sensitivity = mean_sens) %>%
  mutate(Cohort = "HBC (PE)")

# 绘图 P3
p3 <- ggplot(data_HBC, aes(x = fpr, y = sensitivity, color = ft_type)) +
  facet_wrap(~Cohort) +
  geom_line(linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = type_col) +
  # [修改 1] y = NULL 去掉标题
  labs(x = "1 - Specificity", y = NULL) +
  common_theme +
  theme(
    strip.background = element_rect(fill = "#63599b", color = NA),
    # [修改 2] 去掉 Y 轴刻度文字和刻度线
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )

# ==============================================================================
# 5. 拼图并保存
# ==============================================================================

combined_plot <- p1 + p2 + p3

print(combined_plot)

# 保存图片
# 注意：因为三张图拼在一起，且去掉了Y轴，可能需要把宽度设大一点，以免挤压
# 我把宽度从 5 调整到了 12，您可以根据需要改回
ggsave('Anti_graph/Fig5_prediction_validation_PE.pdf', combined_plot, width = 10, height = 4.2)





rm(list = ls())
setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code') 

library(dplyr)
library(ggplot2)
library(readxl)
library(stringr)
library(patchwork)

# ==============================================================================
# 1. 全局设置
# ==============================================================================

N_repeats <- 10 
target_models <- c("Baseline", "Met+Fun", "Met+Bac", "Multi-Omics")

type_col <- c(
  "Baseline" = "#bebebe",
  "Met+Fun" = "#B88D8E",
  "Met+Bac" = "#F0BD52",
  "Multi-Omics" = "#CF7F7D"
)

fix_names <- function(x) {
  x <- gsub("\\.", "+", x) 
  x <- gsub("Multi\\+Omics", "Multi-Omics", x) 
  return(x)
}

calc_signif_label <- function(mean_treat, sd_treat, mean_ctrl, sd_ctrl, n = N_repeats) {
  se_treat <- sd_treat / sqrt(n)
  se_ctrl <- sd_ctrl / sqrt(n)
  se_diff <- sqrt(se_treat^2 + se_ctrl^2)
  t_stat <- (mean_treat - mean_ctrl) / se_diff
  df <- 2 * n - 2
  p_val <- 2 * pt(-abs(t_stat), df)
  
  if (is.na(p_val)) return("")
  if (p_val < 0.001) return("***")
  if (p_val < 0.01) return("**")
  if (p_val < 0.05) return("*")
  return("ns")
}

# ==============================================================================
# 2. 数据处理 (保持不变)
# ==============================================================================

df_thsbc <- read_excel('Anti_results/prediction_summary_auc_THSBC.xlsx') %>%
  mutate(ft_type = fix_names(ft_type)) %>% filter(outcome %in% c('PE_vs_NP', 'PE'), ft_type %in% target_models) %>%
  select(ft_type, roc_auc, auc_sd) %>% rename(mean_auc = roc_auc, sd_auc = auc_sd) %>% mutate(Cohort = "THSBC")

df_webirth <- read_excel('Anti_results/Final_Summary_Stats_WeBirth.xlsx') %>%
  mutate(ft_type = fix_names(ft_type)) %>% filter(outcome %in% c('PE_vs_NP', 'PE'), ft_type %in% target_models) %>%
  group_by(ft_type) %>% slice_max(mean_auc, n = 1) %>% ungroup() %>%
  select(ft_type, mean_auc, sd_auc) %>% mutate(Cohort = "WeBirth")

df_hbc <- read_excel('Anti_results/Final_Summary_Stats_HBC.xlsx') %>%
  mutate(ft_type = fix_names(ft_type)) %>% filter(outcome %in% c('PE_vs_NP', 'PE'), ft_type %in% target_models) %>%
  group_by(ft_type) %>% slice_max(mean_auc, n = 1) %>% ungroup() %>%
  select(ft_type, mean_auc, sd_auc) %>% mutate(Cohort = "HBC")

plot_data <- bind_rows(df_thsbc, df_webirth, df_hbc) %>%
  mutate(ft_type = factor(ft_type, levels = target_models)) %>%
  mutate(Cohort = factor(Cohort, levels = c("THSBC", "WeBirth", "HBC"))) %>%
  group_by(Cohort) %>%
  mutate(
    base_mean = mean_auc[ft_type == "Baseline"],
    base_sd = sd_auc[ft_type == "Baseline"],
    signif_label = mapply(calc_signif_label, mean_auc, sd_auc, base_mean, base_sd),
    signif_label = ifelse(ft_type == "Baseline", "", signif_label)
  ) %>%
  ungroup()

# ==============================================================================
# 3. 绘图 (修改了 geom_text 的位置)
# ==============================================================================

p <- ggplot(plot_data, aes(x = Cohort, y = mean_auc, fill = ft_type)) +
  
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.7) +
  
  geom_errorbar(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc),
                width = 0.2, position = position_dodge(width = 0.8), linewidth = 0.4) +
  
  # [修改点 1] 均值数值：让 y 轴位置等于 mean_auc + sd_auc
  geom_text(aes(y = mean_auc + sd_auc, label = sprintf("%.2f", mean_auc)), 
            position = position_dodge(width = 0.8), 
            vjust = -0.5, # 在误差棒上方一点点
            size = 4.3, fontface = "bold") +
  
  # [修改点 2] 显著性标签：也基于 mean_auc + sd_auc，但推得更高
  geom_text(aes(y = mean_auc + sd_auc, label = signif_label), 
            position = position_dodge(width = 0.8), 
            vjust = -1.5, # 推得比数值更高，避免重叠
            size = 5, color = "black", fontface = "bold") +
  
  scale_fill_manual(values = type_col) +
  
  # [修改点 3] 稍微调高 y 轴上限，防止星星被切掉
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  
  labs(y = "AUC (Mean ± SD)", x = NULL, fill = NULL) +
  
  theme_bw() +
  theme(
    panel.border = element_blank(),
    
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.y = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.text.x = element_text(size = 22, color = "black"),
    axis.ticks.x = element_blank(),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    legend.position = c(0.85, 0.92),
    # legend.justification = c(1, 1),
    legend.background = element_blank(), # 去除图例背景色(透明)
    legend.key = element_blank(),        # 去除图例图标的背景灰框
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.5, "cm"),
    
    plot.margin = margin(10, 10, 10, 10)
  )

print(p)

ggsave('Anti_graph/Fig5_Comparison_Barplot_NumbersOnTop.pdf', p, width = 6.5, height = 4.5)









rm(list = ls())
setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code') # 请确保路径正确

library(dplyr)
library(ggplot2)
library(readxl)
library(patchwork)
library(stringr)

# ==============================================================================
# 1. 全局设置：颜色与模型顺序
# ==============================================================================

# 指定模型顺序
target_models <- c("Baseline", "Met+Fun", "Met+Bac", "Multi-Omics")

# 模型线条颜色
type_col <- c(
  "Baseline" = "#bebebe",
  "Met+Fun" = "#9B6B6C",
  "Met+Bac" = "#E69F00",
  "Multi-Omics" = "#B24745"
)

# 定义通用的分面文字样式
strip_text_style <- element_text(size = 16, color = 'white', face = "bold")

# 定义通用主题 (去除背景网格，保留边框)
common_theme <- theme_bw() +
  theme(
    axis.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    strip.text = strip_text_style 
  )

# ==============================================================================
# 2. 处理 THSBC 数据 (左图 - 保留 Y 轴)
# ==============================================================================

data_thsbc <- read_excel('Anti_results/prediction_roc_curve_data_THSBC.xlsx')

plot_data_thsbc <- data_thsbc %>%
  filter(outcome %in% c('Combined_vs_NP', 'PE/GH')) %>%
  filter(ft_type %in% target_models) %>%
  mutate(ft_type = factor(ft_type, levels = target_models)) %>%
  rename(sensitivity = mean_sens) %>%
  mutate(Cohort = "THSBC (PE/GH)")

# 绘图 P1 (保留完整 Y 轴)
p1 <- ggplot(plot_data_thsbc, aes(x = fpr, y = sensitivity, color = ft_type)) +
  facet_wrap(~Cohort) + 
  geom_line(linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = type_col) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  common_theme +
  theme(strip.background = element_rect(fill = "#d66a43", color = NA))

# ==============================================================================
# 3. 处理 WeBirth 数据 (中图 - 去除 Y 轴)
# ==============================================================================

data_webirth <- read.csv('Anti_results/Final_ROC_Curves_Data_WeBirth_combine.csv')

if("sens" %in% colnames(data_webirth)) {
  data_webirth <- data_webirth %>% rename(mean_sens = sens)
}

plot_data_webirth <- data_webirth %>%
  filter(outcome %in% c('Combined_vs_NP', 'PE/GH')) %>%
  filter(ft_type %in% target_models) %>%
  mutate(ft_type = factor(ft_type, levels = target_models)) %>%
  rename(sensitivity = mean_sens) %>%
  mutate(Cohort = "WeBirth (PE/GH)")

# 绘图 P2
p2 <- ggplot(plot_data_webirth, aes(x = fpr, y = sensitivity, color = ft_type)) +
  facet_wrap(~Cohort) +
  geom_line(linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = type_col) +
  # [修改 1] y = NULL 去掉标题
  labs(x = "1 - Specificity", y = NULL) + 
  common_theme +
  theme(
    strip.background = element_rect(fill = "#a6bfbd", color = NA),
    # [修改 2] 去掉 Y 轴刻度文字和刻度线
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )

# ==============================================================================
# 4. 处理 HBC 数据 (右图 - 去除 Y 轴)
# ==============================================================================

data_HBC <- read.csv('Anti_results/Final_ROC_Curves_Data_HBC.csv')

if("sens" %in% colnames(data_HBC)) {
  data_HBC <- data_HBC %>% rename(mean_sens = sens)
}

data_HBC <- data_HBC %>%
  filter(outcome %in% c('Combined_vs_NP', 'PE/GH')) %>%
  filter(ft_type %in% target_models) %>%
  mutate(ft_type = factor(ft_type, levels = target_models)) %>%
  rename(sensitivity = mean_sens) %>%
  mutate(Cohort = "HBC (PE/GH)")

# 绘图 P3
p3 <- ggplot(data_HBC, aes(x = fpr, y = sensitivity, color = ft_type)) +
  facet_wrap(~Cohort) +
  geom_line(linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = type_col) +
  # [修改 1] y = NULL 去掉标题
  labs(x = "1 - Specificity", y = NULL) +
  common_theme +
  theme(
    strip.background = element_rect(fill = "#63599b", color = NA),
    # [修改 2] 去掉 Y 轴刻度文字和刻度线
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )

# ==============================================================================
# 5. 拼图并保存
# ==============================================================================

combined_plot <- p1 + p2 + p3

print(combined_plot)

# 保存图片
# 注意：因为三张图拼在一起，且去掉了Y轴，可能需要把宽度设大一点，以免挤压
# 我把宽度从 5 调整到了 12，您可以根据需要改回
ggsave('Anti_graph/Figs5_prediction_validation_combine.pdf', combined_plot, width = 10, height = 4.2)





rm(list = ls())
setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code') 

library(dplyr)
library(ggplot2)
library(readxl)
library(stringr)
library(patchwork)

# ==============================================================================
# 1. 全局设置
# ==============================================================================

N_repeats <- 10 
target_models <- c("Baseline", "Met+Fun", "Met+Bac", "Multi-Omics")

type_col <- c(
  "Baseline" = "#bebebe",
  "Met+Fun" = "#B88D8E",
  "Met+Bac" = "#F0BD52",
  "Multi-Omics" = "#CF7F7D"
)

fix_names <- function(x) {
  x <- gsub("\\.", "+", x) 
  x <- gsub("Multi\\+Omics", "Multi-Omics", x) 
  return(x)
}

calc_signif_label <- function(mean_treat, sd_treat, mean_ctrl, sd_ctrl, n = N_repeats) {
  se_treat <- sd_treat / sqrt(n)
  se_ctrl <- sd_ctrl / sqrt(n)
  se_diff <- sqrt(se_treat^2 + se_ctrl^2)
  t_stat <- (mean_treat - mean_ctrl) / se_diff
  df <- 2 * n - 2
  p_val <- 2 * pt(-abs(t_stat), df)
  
  if (is.na(p_val)) return("")
  if (p_val < 0.001) return("***")
  if (p_val < 0.01) return("**")
  if (p_val < 0.05) return("*")
  return("ns")
}

# ==============================================================================
# 2. 数据处理 (保持不变)
# ==============================================================================

df_thsbc <- read_excel('Anti_results/prediction_summary_auc_THSBC.xlsx') %>%
  mutate(ft_type = fix_names(ft_type)) %>% filter(outcome %in% c('Combined_vs_NP', 'PE'), ft_type %in% target_models) %>%
  select(ft_type, roc_auc, auc_sd) %>% rename(mean_auc = roc_auc, sd_auc = auc_sd) %>% mutate(Cohort = "THSBC")

df_webirth <- read_excel('Anti_results/Final_Summary_Stats_WeBirth_combine.xlsx') %>%
  mutate(ft_type = fix_names(ft_type)) %>% filter(outcome %in% c('Combined_vs_NP', 'PE'), ft_type %in% target_models) %>%
  group_by(ft_type) %>% slice_max(mean_auc, n = 1) %>% ungroup() %>%
  select(ft_type, mean_auc, sd_auc) %>% mutate(Cohort = "WeBirth")

df_hbc <- read_excel('Anti_results/Final_Summary_Stats_HBC.xlsx') %>%
  mutate(ft_type = fix_names(ft_type)) %>% filter(outcome %in% c('Combined_vs_NP', 'PE'), ft_type %in% target_models) %>%
  group_by(ft_type) %>% slice_max(mean_auc, n = 1) %>% ungroup() %>%
  select(ft_type, mean_auc, sd_auc) %>% mutate(Cohort = "HBC")

plot_data <- bind_rows(df_thsbc, df_webirth, df_hbc) %>%
  mutate(ft_type = factor(ft_type, levels = target_models)) %>%
  mutate(Cohort = factor(Cohort, levels = c("THSBC", "WeBirth", "HBC"))) %>%
  group_by(Cohort) %>%
  mutate(
    base_mean = mean_auc[ft_type == "Baseline"],
    base_sd = sd_auc[ft_type == "Baseline"],
    signif_label = mapply(calc_signif_label, mean_auc, sd_auc, base_mean, base_sd),
    signif_label = ifelse(ft_type == "Baseline", "", signif_label)
  ) %>%
  ungroup()

# ==============================================================================
# 3. 绘图 (修改了 geom_text 的位置)
# ==============================================================================

p <- ggplot(plot_data, aes(x = Cohort, y = mean_auc, fill = ft_type)) +
  
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.7) +
  
  geom_errorbar(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc),
                width = 0.2, position = position_dodge(width = 0.8), linewidth = 0.4) +
  
  # [修改点 1] 均值数值：让 y 轴位置等于 mean_auc + sd_auc
  geom_text(aes(y = mean_auc + sd_auc, label = sprintf("%.2f", mean_auc)), 
            position = position_dodge(width = 0.8), 
            vjust = -0.5, # 在误差棒上方一点点
            size = 4.3, fontface = "bold") +
  
  # [修改点 2] 显著性标签：也基于 mean_auc + sd_auc，但推得更高
  geom_text(aes(y = mean_auc + sd_auc, label = signif_label), 
            position = position_dodge(width = 0.8), 
            vjust = -1.5, # 推得比数值更高，避免重叠
            size = 5, color = "black", fontface = "bold") +
  
  scale_fill_manual(values = type_col) +
  
  # [修改点 3] 稍微调高 y 轴上限，防止星星被切掉
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  
  labs(y = "AUC (Mean ± SD)", x = NULL, fill = NULL) +
  
  theme_bw() +
  theme(
    panel.border = element_blank(),
    
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.y = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.text.x = element_text(size = 22, color = "black"),
    axis.ticks.x = element_blank(),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    legend.position = c(0.85, 0.92),
    # legend.justification = c(1, 1),
    legend.background = element_blank(), # 去除图例背景色(透明)
    legend.key = element_blank(),        # 去除图例图标的背景灰框
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.5, "cm"),
    
    plot.margin = ggplot2::margin(10, 10, 10, 10)
  )

print(p)

ggsave('Anti_graph/Figs5_Comparison_Barplot_NumbersOnTop.pdf', p, width = 6.5, height = 4.5)


