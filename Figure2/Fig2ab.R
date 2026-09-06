rm(list = ls())
setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code')

library(dplyr)
library(data.table)
library(stringr)
library(reshape2)
library(openxlsx)
library(ggplot2)
library(ggh4x)
library(vegan)
library(ggpubr)
library(patchwork)

library(lme4); library(lmerTest); library(pbkrtest)

#############################Richness################################




data = read.csv('data/Weighted_ALL_ITS.csv')

data <- data[!is.na(data$wt1), ]

table(data$preeclampsia)

data$period <- sub("^V", "T", as.character(data$period))

data <- data %>%
  rename(PE = preeclampsia) %>%
  filter(!is.na(HDP)) %>%
  mutate(
    GH = ifelse(HDP == 1 & PE == 0, 1, 0),
    index = "Richness"
  )
taxa_cols <- grep("^g.*[0-9]+$", names(data), value = TRUE)
data$alpha <- rowSums(data[, taxa_cols] > 0, na.rm = TRUE)

datac <- data %>%
  dplyr::select(id, period, HDP, PE, GH, index, alpha)

table(datac$period)

length(unique(datac$id))

data_1 = read.csv('data/Weighted_ALL_Taxonomy.csv')

data_1 <- data_1[!is.na(data_1$wt2), ]


data_1$period <- sub("^V", "T", as.character(data_1$period))

data_1 <- data_1 %>%
  rename(PE = preeclampsia) %>%
  filter(!is.na(HDP)) %>%
  mutate(
    GH = ifelse(HDP == 1 & PE == 0, 1, 0),
    index = "Richness"
  )
taxa_cols <- grep("^s.*[0-9]+$", names(data_1), value = TRUE)
data_1$alpha <- rowSums(data_1[, taxa_cols] > 0, na.rm = TRUE)

datak <- data_1 %>%
  select(id, period, HDP, PE, GH, index, alpha)




datac$source <- "ITS" #


datak$source <- "Metagenomic species" #




data_16S <-  datac
data_16S$source <- "ITS" #

data_16S <- data_16S %>%
  mutate(
    group = case_when(
      PE == 1 ~ "PE",
      GH == 1 ~ "GH",
      TRUE    ~ "NP"
    )
  )



filtered_data <- data_16S[data_16S$index == "Richness", ]


names(filtered_data)[names(filtered_data) == "alpha"] <- "Shannon"



# calculate_trend_pv <- function(df, group_var = "group", metric = "Shannon", ref = "NP") {
#   df <- df %>%
#     mutate(
#       period = as.numeric(str_remove(period, "T")),
#       !!group_var := factor(.data[[group_var]], levels = c("NP", "PE", "GH"))
#     )
#   
#   # 每组单独模型：metric ~ (1|id) + period
#   for (g in levels(df[[group_var]])) {
#     fit_g <- lmer(as.formula(paste0(metric, " ~ (1|id) + period")),
#                   data = df[df[[group_var]] == g, ], REML = FALSE)
#     cat("\nGroup:", g, "\n")
#     print(round(summary(fit_g)$coefficients["period", ], 4))
#   }
#   
#   # 整体交互模型：metric ~ (1|id) + group * period
#   df[[group_var]] <- relevel(df[[group_var]], ref = ref)
#   fit_all <- lmer(as.formula(paste0(metric, " ~ (1|id) + ", group_var, " * period")),
#                   data = df, REML = FALSE)
#   coefs <- summary(fit_all)$coefficients
#   
#   cat("\nOverall model (reference =", ref, "):\n")
#   # 参考组的 period 效应
#   if ("period" %in% rownames(coefs)) {
#     cat("Reference period effect:\n")
#     print(round(coefs["period", ], 4))
#   }
#   # 交互项（各组与参考组的斜率差）
#   for (g in setdiff(levels(df[[group_var]]), ref)) {
#     term <- paste0(group_var, g, ":period")
#     if (term %in% rownames(coefs)) {
#       cat("Interaction (", g, " vs ", ref, "):\n", sep = "")
#       print(round(coefs[term, ], 4))
#     }
#   }
#   
#   # 新增：计算并输出 “Interaction (PE vs GH)”
#   # 通过把参考组改为 GH，读取 groupPE:period 即得到 PE - GH 的斜率差
#   df_gh <- df
#   df_gh[[group_var]] <- relevel(df_gh[[group_var]], ref = "GH")
#   fit_pe_vs_gh <- lmer(as.formula(paste0(metric, " ~ (1|id) + ", group_var, " * period")),
#                        data = df_gh, REML = FALSE)
#   coefs_gh <- summary(fit_pe_vs_gh)$coefficients
#   term_pe_vs_gh <- paste0(group_var, "PE:period")
#   if (term_pe_vs_gh %in% rownames(coefs_gh)) {
#     cat("\nInteraction (PE vs GH):\n")
#     print(round(coefs_gh[term_pe_vs_gh, ], 4))
#   }
#   
#   invisible(list(
#     fit_all_ref = fit_all,
#     fit_pe_vs_gh = fit_pe_vs_gh
#   ))
# }

library(dplyr)
library(stringr)


calculate_trend_pv <- function(df, group_var = "group", metric = "Shannon", ref = "NP") {
  df <- df %>%
    mutate(
      period = as.numeric(str_remove(period, "T")),
      !!group_var := factor(.data[[group_var]], levels = c("NP", "PE", "GH"))
    )
  
  # 创建一个空的数据框，用于收集后续所有的统计结果
  res_list <- list()
  
  # 1. 每组单独模型：metric ~ (1|id) + period
  for (g in levels(df[[group_var]])) {
    fit_g <- lmer(as.formula(paste0(metric, " ~ (1|id) + period")),
                  data = df[df[[group_var]] == g, ], REML = FALSE)
    
    coefs <- summary(fit_g)$coefficients
    # lmerTest 提供的系数表最后一列通常是 p 值 (Pr(>|t|))
    p_val <- coefs["period", ncol(coefs)] 
    
    res_list[[paste0(g, "_trend")]] <- data.frame(
      Test = paste0(g, " internal trend"),
      Estimate = coefs["period", "Estimate"],
      P_value = p_val
    )
  }
  
  # 2. 整体交互模型：metric ~ (1|id) + group * period
  df[[group_var]] <- relevel(df[[group_var]], ref = ref)
  fit_all <- lmer(as.formula(paste0(metric, " ~ (1|id) + ", group_var, " * period")),
                  data = df, REML = FALSE)
  coefs_all <- summary(fit_all)$coefficients
  
  # 提取交互项（各组与参考组的斜率差）
  for (g in setdiff(levels(df[[group_var]]), ref)) {
    term <- paste0(group_var, g, ":period")
    if (term %in% rownames(coefs_all)) {
      res_list[[paste0(g, "_vs_", ref)]] <- data.frame(
        Test = paste0("Interaction (", g, " vs ", ref, ")"),
        Estimate = coefs_all[term, "Estimate"],
        P_value = coefs_all[term, ncol(coefs_all)]
      )
    }
  }
  
  # 3. 计算并提取 “Interaction (PE vs GH)”
  df_gh <- df
  df_gh[[group_var]] <- relevel(df_gh[[group_var]], ref = "GH")
  fit_pe_vs_gh <- lmer(as.formula(paste0(metric, " ~ (1|id) + ", group_var, " * period")),
                       data = df_gh, REML = FALSE)
  coefs_gh <- summary(fit_pe_vs_gh)$coefficients
  term_pe_vs_gh <- paste0(group_var, "PE:period")
  
  if (term_pe_vs_gh %in% rownames(coefs_gh)) {
    res_list[["PE_vs_GH"]] <- data.frame(
      Test = "Interaction (PE vs GH)",
      Estimate = coefs_gh[term_pe_vs_gh, "Estimate"],
      P_value = coefs_gh[term_pe_vs_gh, ncol(coefs_gh)]
    )
  }
  
  # 4. 合并所有结果并进行 FDR 矫正 -------------------
  final_results <- bind_rows(res_list)
  # 使用自带的 p.adjust 函数进行 FDR 矫正 (Benjamini-Hochberg 方法)
  final_results$FDR <- p.adjust(final_results$P_value, method = "fdr")
  
  # 打印漂亮的汇总表格
  cat("\n=== Trend & Interaction Statistics ===\n")
  print(final_results %>% mutate_if(is.numeric, round, 4)) # 保留4位小数输出
  cat("======================================\n")
  
  # 返回模型和统计表
  invisible(list(
    stats_table = final_results,
    fit_all_ref = fit_all,
    fit_pe_vs_gh = fit_pe_vs_gh
  ))
}


calculate_trend_pv(filtered_data, group_var = "group", metric = "Shannon", ref = "NP")






col_group = c('PE'='#cc8082', "GH" = "#f9ddba", 'NP'='#5c6d8e')

filtered_data <- filtered_data %>%
  mutate(group = factor(group, levels = c("NP", "GH", "PE")))

my_comparisons <- list(
  c("NP", "GH"),
  c("NP", "PE"),
  c("GH", "PE")
)




# p_shannon <- ggplot(filtered_data, aes(x = group, y = Shannon, fill = group, color = group)) +
#   facet_wrap(~ period, nrow = 1) +
#   geom_violin(width = 0.3, alpha = 0.2, position = position_dodge(0.4), trim = TRUE) +
#   geom_boxplot(size = 0.8, width = 0.3, alpha = 0.1, outlier.shape = NA,
#                position = position_dodge(0.4)) +
#   stat_summary(fun = median, geom = "line", aes(group = group),
#                position = position_dodge(0.4)) +
#   scale_fill_manual(values = col_group) +
#   scale_color_manual(values = col_group) +
#   scale_y_continuous(expand = c(0, 0.1), limits = c(NA, 105)) +
#   theme_bw() +
#   theme(
#     axis.title = element_blank(),
#     axis.text.x = element_text(size = 11.2, color = "black"),
#     axis.text.y = element_text(size = 13.2, color = "grey10"),
#     strip.text = element_text(size = 11, color = 'black'),
#     strip.background = element_rect(fill = '#c0a8b6'),
#     legend.title = element_blank(),
#     panel.grid = element_line(linetype = 2),
#     legend.position = "none"
#   ) +
#   stat_compare_means(
#     comparisons = my_comparisons,
#     method = "wilcox.test",      # 或 "t.test"，看你数据分布
#     label = "p.signif",          # 或 "p.format"
#     hide.ns = FALSE              # TRUE 可隐藏不显著结果
#   )


library(rstatix)
library(ggpubr)

# 第一步：在画图前，先单独计算并进行 FDR 矫正

stat.test <- filtered_data %>%   
  group_by(period) %>%   
  wilcox_test(Shannon ~ group, comparisons = my_comparisons) %>%   
  group_by(period) %>%   
  adjust_pvalue(method = "fdr") %>%   
  add_significance("p.adj") %>%   
  ungroup()

stat.test <- stat.test %>%
  add_xy_position(
    x = "group",
    fun = "max",
    data = filtered_data,
    formula = Shannon ~ group
  )

# 第二步：原来的画图代码，只替换最后一句
p_shannon <- ggplot(filtered_data, aes(x = group, y = Shannon, fill = group, color = group)) +
  facet_wrap(~ period, nrow = 1) +
  geom_violin(width = 0.3, alpha = 0.2, position = position_dodge(0.4), trim = TRUE) +
  geom_boxplot(size = 0.8, width = 0.3, alpha = 0.1, outlier.shape = NA,
               position = position_dodge(0.4)) +
  stat_summary(fun = median, geom = "line", aes(group = group),
               position = position_dodge(0.4)) +
  scale_fill_manual(values = col_group) +
  scale_color_manual(values = col_group) +
  scale_y_continuous(expand = c(0, 0.1), limits = c(NA, 300)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 11.2, color = "black"),
    axis.text.y = element_text(size = 13.2, color = "grey10"),
    strip.text = element_text(size = 11, color = 'black'),
    strip.background = element_rect(fill = '#c0a8b6'),
    legend.title = element_blank(),
    panel.grid = element_line(linetype = 2),
    legend.position = "none"
  ) +
  stat_pvalue_manual(
    stat.test,               
    label = "p.adj",  ###label = "p.adj.signif"
    tip.length = 0.01,
    hide.ns = FALSE,
    inherit.aes = FALSE      # <--- 加上这一行，断开全局继承，直接解决报错！
  )

print(stat.test)



p_shannon <- ggplot(filtered_data, aes(x = period, y = Shannon, fill = group, color = group)) +
  facet_wrap(~ "Fungal Richness") +
  geom_violin(width = 0.3, alpha = 0.2, position = position_dodge(0.4), trim = TRUE) +
  geom_boxplot(size = 0.8, width = 0.3, alpha = 0.1, outlier.shape = NA,
               position = position_dodge(0.4)) +
  stat_summary(fun = median, geom = "line", aes(group = group),
               position = position_dodge(0.4)) +
  scale_fill_manual(values = col_group) +
  scale_color_manual(values = col_group) +
  scale_y_continuous(expand = c(0, 0.1), limits = c(NA, 80)) +
  annotate('text', x = 1.6, y = 78, size = 4., hjust = 0,
           label=expression(paste(italic('P'),'-interaction (PE vs NP) = 0.067')))+
  annotate('text', x = 1.6, y = 75, size = 4., hjust = 0,
           label=expression(paste(italic('P'),'-interaction (GH vs NP) < 0.001')))+
  annotate('text', x = 1.6, y = 72, size = 4., hjust = 0,
           label=expression(paste(italic('P'),'-interaction (PE vs GH) = 0.052')))+
  annotate('text', x = 1.6, y = 69, size = 4., hjust = 0,
           label=expression(paste(italic('P'),'-trend < 0.001 for PE')))+
  annotate('text', x = 1.6, y = 66, size = 4., hjust = 0,
           label=expression(paste(italic('P'),'-trend < 0.001 for GH')))+
  annotate('text', x = 1.6, y = 63, size = 4., hjust = 0,
           label=expression(paste(italic('P'),'-trend < 0.001 for NP')))+
  geom_signif(
    # 每个period的三个比较，按照新顺序
    xmin = c(1.0, 0.85, 0.85,     # period 1: GH vs Non-HDP, PE vs GH, PE vs Non-HDP
             2.0, 1.85, 1.85,      # period 2: 同样顺序
             3.0, 2.85, 2.85),     # period 3: 同样顺序
    xmax = c(1.15, 1.0, 1.15,      
             2.15, 2.0, 2.15,      
             3.15, 3.0, 3.15),  
    y_position = c(62, 66, 70,  # period 1的不同高度
                   47, 51, 55,  # period 2的不同高度
                   47, 51, 55), # period 3的不同高度
    annotations = c("NS", "***", "NS",    # period 1的三个比较
                    "NS", "NS", "NS",    # period 2的三个比较
                    "NS", "NS", "NS"), # 根据您的实际显著性结果调整
    tip_length = 0.01, 
    size = 0.4, 
    textsize = 3.7, 
    vjust = 0.05
  )+
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 17, color = "black"),
    axis.text.y = element_text(size = 18, color = "grey10"),
    strip.text = element_text(size = 20, color = 'black'),
    strip.background = element_rect(fill = '#e6e4c5'),
    legend.title = element_blank(),
    panel.grid = element_line(linetype = 2),
    legend.position = "none"
  )



ggsave('Anti_graph/fig2a_Richness_ITS_class3.pdf', p_shannon, width = 4.6, height = 4.3)






###################Richness species############################


data_16S <-  datak
range(datak$alpha, na.rm = TRUE)
data_16S$source <- "16S"

# data_16S$group <- ifelse(data_16S$HDP == 1, "HDP", "Non-HDP")


data_16S <- data_16S %>%
  mutate(
    group = case_when(
      PE == 1 ~ "PE",
      GH == 1 ~ "GH",
      TRUE    ~ "NP"
    )
  )

filtered_data <- data_16S[data_16S$index == "Richness", ]

names(filtered_data)[names(filtered_data) == "alpha"] <- "Shannon"




# calculate_trend_pv <- function(df, group_var = "group", metric = "Shannon", ref = "NP") {
#   df <- df %>%
#     mutate(
#       period = as.numeric(str_remove(period, "T")),
#       !!group_var := factor(.data[[group_var]], levels = c("NP", "PE", "GH"))
#     )
#   
#   # 每组单独模型：metric ~ (1|id) + period
#   for (g in levels(df[[group_var]])) {
#     fit_g <- lmer(as.formula(paste0(metric, " ~ (1|id) + period")),
#                   data = df[df[[group_var]] == g, ], REML = FALSE)
#     cat("\nGroup:", g, "\n")
#     print(round(summary(fit_g)$coefficients["period", ], 4))
#   }
#   
#   # 整体交互模型：metric ~ (1|id) + group * period
#   df[[group_var]] <- relevel(df[[group_var]], ref = ref)
#   fit_all <- lmer(as.formula(paste0(metric, " ~ (1|id) + ", group_var, " * period")),
#                   data = df, REML = FALSE)
#   coefs <- summary(fit_all)$coefficients
#   
#   cat("\nOverall model (reference =", ref, "):\n")
#   # 参考组的 period 效应
#   if ("period" %in% rownames(coefs)) {
#     cat("Reference period effect:\n")
#     print(round(coefs["period", ], 4))
#   }
#   # 交互项（各组与参考组的斜率差）
#   for (g in setdiff(levels(df[[group_var]]), ref)) {
#     term <- paste0(group_var, g, ":period")
#     if (term %in% rownames(coefs)) {
#       cat("Interaction (", g, " vs ", ref, "):\n", sep = "")
#       print(round(coefs[term, ], 4))
#     }
#   }
#   
#   # 新增：计算并输出 “Interaction (PE vs GH)”
#   # 通过把参考组改为 GH，读取 groupPE:period 即得到 PE - GH 的斜率差
#   df_gh <- df
#   df_gh[[group_var]] <- relevel(df_gh[[group_var]], ref = "GH")
#   fit_pe_vs_gh <- lmer(as.formula(paste0(metric, " ~ (1|id) + ", group_var, " * period")),
#                        data = df_gh, REML = FALSE)
#   coefs_gh <- summary(fit_pe_vs_gh)$coefficients
#   term_pe_vs_gh <- paste0(group_var, "PE:period")
#   if (term_pe_vs_gh %in% rownames(coefs_gh)) {
#     cat("\nInteraction (PE vs GH):\n")
#     print(round(coefs_gh[term_pe_vs_gh, ], 4))
#   }
#   
#   invisible(list(
#     fit_all_ref = fit_all,
#     fit_pe_vs_gh = fit_pe_vs_gh
#   ))
# }



calculate_trend_pv(filtered_data, group_var = "group", metric = "Shannon", ref = "NP")

# calculate_trend_pv(filtered_data, 'HDP', 'Shannon')

# col_group = c('HDP'='#cc8082','Non-HDP'='#5c6d8e')

col_group = c('PE'='#cc8082', "GH" = "#f9ddba", 'NP'='#5c6d8e')




filtered_data <- filtered_data %>%
  mutate(group = factor(group, levels = c("NP", "GH", "PE")))

my_comparisons <- list(
  c("NP", "GH"),
  c("NP", "PE"),
  c("GH", "PE")
)



# p_shannon <- ggplot(filtered_data, aes(x = group, y = Shannon, fill = group, color = group)) +
#   facet_wrap(~ period, nrow = 1) +
#   geom_violin(width = 0.3, alpha = 0.2, position = position_dodge(0.4), trim = TRUE) +
#   geom_boxplot(size = 0.8, width = 0.3, alpha = 0.1, outlier.shape = NA,
#                position = position_dodge(0.4)) +
#   stat_summary(fun = median, geom = "line", aes(group = group),
#                position = position_dodge(0.4)) +
#   scale_fill_manual(values = col_group) +
#   scale_color_manual(values = col_group) +
#   scale_y_continuous(expand = c(0, 0.1), limits = c(NA, 355)) +
#   theme_bw() +
#   theme(
#     axis.title = element_blank(),
#     axis.text.x = element_text(size = 11.2, color = "black"),
#     axis.text.y = element_text(size = 13.2, color = "grey10"),
#     strip.text = element_text(size = 11, color = 'black'),
#     strip.background = element_rect(fill = '#c0a8b6'),
#     legend.title = element_blank(),
#     panel.grid = element_line(linetype = 2),
#     legend.position = "none"
#   ) +
#   stat_compare_means(
#     comparisons = my_comparisons,
#     method = "wilcox.test",      # 或 "t.test"，看你数据分布
#     label = "p.signif",          # 或 "p.format"
#     hide.ns = FALSE              # TRUE 可隐藏不显著结果
#   )


library(rstatix)
library(ggpubr)

# 第一步：在画图前，先单独计算并进行 FDR 矫正
stat.test <- filtered_data %>%   
  group_by(period) %>%   
  wilcox_test(Shannon ~ group, comparisons = my_comparisons) %>%   
  group_by(period) %>%   
  adjust_pvalue(method = "fdr") %>%   
  add_significance("p.adj") %>%   
  ungroup()

stat.test <- stat.test %>%
  add_xy_position(
    x = "group",
    fun = "max",
    data = filtered_data,
    formula = Shannon ~ group
  )

# 第二步：原来的画图代码，只替换最后一句
p_shannon <- ggplot(filtered_data, aes(x = group, y = Shannon, fill = group, color = group)) +
  facet_wrap(~ period, nrow = 1) +
  geom_violin(width = 0.3, alpha = 0.2, position = position_dodge(0.4), trim = TRUE) +
  geom_boxplot(size = 0.8, width = 0.3, alpha = 0.1, outlier.shape = NA,
               position = position_dodge(0.4)) +
  stat_summary(fun = median, geom = "line", aes(group = group),
               position = position_dodge(0.4)) +
  scale_fill_manual(values = col_group) +
  scale_color_manual(values = col_group) +
  scale_y_continuous(expand = c(0, 0.1), limits = c(NA, 300)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 11.2, color = "black"),
    axis.text.y = element_text(size = 13.2, color = "grey10"),
    strip.text = element_text(size = 11, color = 'black'),
    strip.background = element_rect(fill = '#c0a8b6'),
    legend.title = element_blank(),
    panel.grid = element_line(linetype = 2),
    legend.position = "none"
  ) +
  stat_pvalue_manual(
    stat.test,               
    label = "p.adj",  ###label = "p.adj.signif"
    tip.length = 0.01,
    hide.ns = FALSE,
    inherit.aes = FALSE      # <--- 加上这一行，断开全局继承，直接解决报错！
  )

print(stat.test)


p_shannon <- ggplot(filtered_data, aes(x = period, y = Shannon, fill = group, color = group)) +
  facet_wrap(~ "Bacterial Richness") +
  geom_violin(width = 0.3, alpha = 0.2, position = position_dodge(0.4), trim = TRUE) +
  geom_boxplot(size = 0.8, width = 0.3, alpha = 0.1, outlier.shape = NA,
               position = position_dodge(0.4)) +
  stat_summary(fun = median, geom = "line", aes(group = group),
               position = position_dodge(0.4)) +
  scale_fill_manual(values = col_group) +
  scale_color_manual(values = col_group) +
  scale_y_continuous(expand = c(0, 0.1), limits = c(NA, 320)) +
  annotate('text', x = 1.5, y = 313, size = 4., hjust = 0,
           label=expression(paste(italic('P'),'-interaction (PE vs NP) < 0.001')))+
  annotate('text', x = 1.5, y = 302, size = 4., hjust = 0,
           label=expression(paste(italic('P'),'-interaction (GH vs NP) = 0.442')))+
  annotate('text', x = 1.5, y = 290, size = 4., hjust = 0,
           label=expression(paste(italic('P'),'-interaction (PE vs GH) = 0.017')))+
  annotate('text', x = 1.5, y = 279, size = 4., hjust = 0,
           label=expression(paste(italic('P'),'-trend < 0.001 for PE')))+
  annotate('text', x = 1.5, y = 268, size = 4., hjust = 0,
           label=expression(paste(italic('P'),'-trend = 0.077 for GH')))+
  annotate('text', x = 1.5, y = 257, size = 4., hjust = 0,
           label=expression(paste(italic('P'),'-trend = 0.010 for NP')))+
  geom_signif(
    # 每个period的三个比较，按照新顺序
    xmin = c(1.0, 0.85, 0.85,     # period 1: GH vs Non-HDP, PE vs GH, PE vs Non-HDP
             2.0, 1.85, 1.85,      # period 2: 同样顺序
             3.0, 2.85, 2.85),     # period 3: 同样顺序
    xmax = c(1.15, 1.0, 1.15,      
             2.15, 2.0, 2.15,      
             3.15, 3.0, 3.15),  
    y_position = c(236, 248, 260,  # period 1的不同高度
                   216, 228, 240,  # period 2的不同高度
                   216, 228, 240), # period 3的不同高度
    annotations = c("NS", "NS", "NS",    # period 1的三个比较
                    "NS", "NS", "NS",    # period 2的三个比较
                    "*", "NS", "**"), # 根据您的实际显著性结果调整PEvsGH xxx PEvsNP
    tip_length = 0.01, 
    size = 0.4, 
    textsize = 3.7
  )+
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 17, color = "black"),
    axis.text.y = element_text(size = 18, color = "grey10"),
    strip.text = element_text(size = 20, color = 'black'),
    strip.background = element_rect(fill = '#fbe5e6'),
    legend.title = element_blank(),
    panel.grid = element_line(linetype = 2),
    legend.position = "none"
  )



# p_shannon = ggplot(filtered_data, aes(x = period, y = Shannon, fill = group, color=group))+
#   facet_wrap(~"Shannon")+
#   geom_violin(width = 0.3, alpha = 0.2, position = position_dodge(0.4))+
#   geom_boxplot(size= 0.8, width = 0.3, alpha = 0.1, outlier.shape = NA, 
#                position = position_dodge(0.4))+
#   # geom_point(alpha=0.5, size=1,
#   #            position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0,
#   #                                            dodge.width = 0.4))+
#   stat_summary(fun = median, geom = "line", aes(group = group), 
#                position = position_dodge(0.4))+
#   scale_fill_manual(values = col_group)+
#   scale_color_manual(values = col_group)+
#   # scale_y_continuous(expand = c(0,0.1))+
#   scale_y_continuous(expand = c(0,0.1), limits = c(1, NA))+
#   annotate('text', x = 0.55, y = 1.86, size = 3.35, hjust = 0,
#            label=expression(paste(italic('P'),'-interaction = 0.224')))+
#   annotate('text', x = 0.55, y = 1.525, size = 3.35, hjust = 0,
#            label=expression(paste(italic('P'),'-trend < 0.001 for HDP')))+
#   annotate('text', x = 0.55, y = 1.20, size = 3.35, hjust = 0,
#            label=expression(paste(italic('P'),'-trend < 0.001 for Non-HDP')))+
#   stat_compare_means(aes(group = group), method = "wilcox.test",
#                      label = "p.signif", show.legend = F, label.y.npc = 0.96)+
#   theme_bw()+
#   theme(
#     axis.title = element_blank(),
#     axis.text.x = element_text(size = 10.2, color = "black"), 
#     axis.text.y = element_text(size = 9.2, color = "grey10"),
#     strip.text = element_text(size = 11, color = 'black'),
#     strip.background = element_rect(fill = '#c0a8b6'),
#     legend.title = element_blank(),
#     panel.grid= element_line(linetype = 2),
#     legend.position = "none" 
#   )


ggsave('Anti_graph/fig2a_Richness_mgs_class3.pdf', p_shannon, width = 4.6, height = 4.3)







rm(list = ls())
setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code')

library(openxlsx)
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)

merged_data <- openxlsx::read.xlsx("Anti_results/mycro_becte_data.xlsx")
table(merged_data$period)
length(unique(merged_data$id))
merged_data$group <- factor(merged_data$group, levels = c("NP", "GH", "PE"))

lmm_inter <- lmer(Shannon_MGS ~ Shannon * group + (1|id), data = merged_data)
summary_inter <- summary(lmm_inter)

merged_data_gh <- merged_data
merged_data_gh$group <- relevel(merged_data_gh$group, ref = "GH")
lmm_inter_gh <- lmer(Shannon_MGS ~ Shannon * group + (1|id), data = merged_data_gh)
summary_inter_gh <- summary(lmm_inter_gh)


m_np <- lmer(Shannon_MGS ~ Shannon + (1|id), data = filter(merged_data, group == "NP"))
m_gh <- lmer(Shannon_MGS ~ Shannon + (1|id), data = filter(merged_data, group == "GH"))
m_pe <- lmer(Shannon_MGS ~ Shannon + (1|id), data = filter(merged_data, group == "PE"))

coef_np <- summary(m_np)$coefficients["Shannon", ]
coef_gh <- summary(m_gh)$coefficients["Shannon", ]
coef_pe <- summary(m_pe)$coefficients["Shannon", ]

np_pe = summary_inter$coefficients["Shannon:groupPE", "Pr(>|t|)"]
np_gh = summary_inter$coefficients["Shannon:groupGH", "Pr(>|t|)"]
gh_pe = summary_inter_gh$coefficients["Shannon:groupPE", "Pr(>|t|)"]


inter_p_raw <- c(
  np_pe = summary_inter$coefficients["Shannon:groupPE", "Pr(>|t|)"],
  np_gh = summary_inter$coefficients["Shannon:groupGH", "Pr(>|t|)"],
  gh_pe = summary_inter_gh$coefficients["Shannon:groupPE", "Pr(>|t|)"]
)
inter_p_adj <- p.adjust(inter_p_raw, method = "fdr")


fmt_inter_p <- function(p, comp) {
  if (p < 0.001) {
    return(sprintf("italic(P)*'-interaction (%s) < 0.001'", comp))
  } else {
    return(sprintf("italic(P)*'-interaction (%s) = %.3f'", comp, p))
  }
}


fmt_group_label <- function(grp, beta, p) {
  if (p < 0.001) {
    return(sprintf("'%s: ' * beta == '%.2f, ' * italic(P) * ' < 0.001'", grp, beta))
  } else {
    return(sprintf("'%s: ' * beta == '%.2f, ' * italic(P) * ' = %.3f'", grp, beta, p))
  }
}


anno_data <- data.frame(
  x = c(76, 76, 76, 76, 76, 76), 
  y = c(290, 275, 260, 245, 230, 215),      
  label = c(
    fmt_inter_p(np_pe, "NP vs PE"),
    fmt_inter_p(np_gh, "NP vs GH"),
    fmt_inter_p(gh_pe, "GH vs PE"),
    fmt_group_label("NP", coef_np["Estimate"], coef_np["Pr(>|t|)"]),
    fmt_group_label("GH", coef_gh["Estimate"], coef_gh["Pr(>|t|)"]),
    fmt_group_label("PE", coef_pe["Estimate"], coef_pe["Pr(>|t|)"])
  ),
  group = factor(
    c("PE", "GH", "NP", "NP", "GH", "PE"),
    levels = c("NP", "GH", "PE")
  )
)


merged_data$pred_MGS <- predict(lmm_inter, re.form = NA)

group_cols = c("#2E5266","#FFC107","#F44336")

p_T1 = ggplot(merged_data, 
              aes(x = Shannon, 
                  y = Shannon_MGS, 
                  color = group, 
                  group = group, 
                  fill = group)) +
  geom_point(alpha = 0.8, size = 1.2) + 
  
  geom_line(aes(y = pred_MGS), linewidth = 1.1) + 
  
  scale_color_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  
  geom_text(data = anno_data, 
            aes(x = x, y = y, label = label, color = group),
            inherit.aes = FALSE, 
            hjust = 1,           # 👈 1 代表右对齐
            size = 4.5, 
            parse = TRUE, 
            show.legend = FALSE) +
  
  labs(
    x = "Fungal richness",
    y = "Bacterial richness",
    color = "Group",
    fill  = "Group"
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = rel(1.5), color = "black"),
    axis.text  = element_text(size = rel(1.5), color = "black"),
    axis.line  = element_line(linewidth = 0.4, color = "black"),
    legend.position = c(0.95, 0.05), 
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "transparent", color = NA)
  )


ggsave('Anti_graph/fig2b_Richness_LMM.pdf', p_T1, width = 4.6, height = 4.3, bg = "transparent")





