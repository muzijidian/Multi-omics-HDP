rm(list = ls())
setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code')
library(dplyr)
library(data.table)
library(stringr)
library(reshape2)
library(openxlsx)
library(ggplot2)
library(vegan)
library(lme4)
library(lmerTest)
library(doParallel)
library(patchwork)
source('help_func.R')
select = dplyr::select


data_info = read.xlsx("../data/metabolites_info_624_THSBC_merge_corrected.xlsx")

data_ITS = read.xlsx("Anti_results/ITS_Explaining_Metabolites_Academic_Standard_adj.xlsx")

data_MGS = read.xlsx("Anti_results/MGS_Explaining_Metabolites_Academic_Standard_adj.xlsx")

info_clean <- data_info %>% 
  select(Index, sugg_cmpd_name, Class.I) %>% 
  distinct(Index, .keep_all = TRUE)

data_ITS_annotated <- data_ITS %>% 
  left_join(info_clean, by = c("metabolite" = "Index")) %>% 
  mutate(source = 'Fungi')

data_MGS_annotated <- data_MGS %>% 
  left_join(info_clean, by = c("metabolite" = "Index")) %>% 
  mutate(source = 'Bacteria')

its_all <- data_ITS_annotated %>% 
  filter(period == 'All')

mgs_all <- data_MGS_annotated %>% 
  filter(period == 'All')

df_all <- rbind(its_all, mgs_all)






df_all$R2 <- df_all$R2 * 100

df_all$sugg_cmpd_name <- str_replace_all(df_all$sugg_cmpd_name, regex(" Acid", ignore_case = TRUE), " acid")

name_mapping <- c(
  'Glycerophospho-N-Arachidonoyl Ethanolamine' = 'GP-NAE',
  'trans-resveratrol-3-O-sulfate' = 'Trans-resveratrol-3-O-sulfate',
  '2,4-diacetamino-2,4,6-triphenoxy-D-mannopyranose' = 'Man2NAc4NAc',
  'Phosphatidylethanolamine lyso alkenyl 16:0' = 'LPE(P-16:0)',
  'Glycine deoxycholic acid'= 'GDCA',
  'Glycochenodeoxycholic acid'= 'GCDCA',
  'Deoxycholic acid'= 'DCA',
  'Taurocholic acid'= 'TCA',
  'Taurochenodesoxycholic acid'= 'TCDCA',
  'Glycoursodeoxycholic acid'= 'GUDCA',
  'Glycohyodeoxycholic acid'= 'GHDCA'
)

df_all$sugg_cmpd_name <- ifelse(
  df_all$sugg_cmpd_name %in% names(name_mapping),
  name_mapping[df_all$sugg_cmpd_name],
  df_all$sugg_cmpd_name
)



col_source <- c("Bacteria" = "#FFB2B9", "Fungi" = "#ebdcb2")

class_colors <- c("#FFFFB3", "#8DD3C7", "#BEBADA", "#FB8072", 
                  "#80B1D3", "#FDB462", "#FCCDE5", "#D9D9D9", 
                  "#BC80BD", "#CCEBC5", "#3498DB", "#FEB29B")
class_order <- c("AlAm", "AAD", "BA", "BSD", "CHD", "FA", 
                 "GP", "HRC", "NUC", "OAD", "SPL", "Others")

names(class_colors) <- class_order

stats_df <- df_all %>%
  group_by(sugg_cmpd_name) %>%
  summarise(total_R2 = sum(R2, na.rm = TRUE)) %>%
  # filter(total_R2 > 15) %>%
  arrange(desc(total_R2)) %>%
  slice_head(n = 30)


df_plot <- df_all %>%
  filter(sugg_cmpd_name %in% stats_df$sugg_cmpd_name)

df_plot$sugg_cmpd_name <- factor(df_plot$sugg_cmpd_name, 
                                 levels = stats_df$sugg_cmpd_name)




p_bar <- ggplot(df_plot, aes(x = sugg_cmpd_name, y = R2, fill = source)) +
  geom_col(width = 0.8, color = NA) + # width 控制柱子宽度
  scale_fill_manual(values = col_source) +
  labs(y = "Variance explained (%)", x = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(), 
    axis.line.x = element_blank(),  
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = c(0.45, 0.85),
    axis.ticks.length.x = unit(0, "mm"),
    plot.margin = ggplot2::margin(t = 5, r = 5, b = 0, l = 5)
  )


df_annotation <- df_plot %>%
  select(sugg_cmpd_name, Class.I) %>%
  distinct()



p_strip <- ggplot(df_annotation, aes(x = sugg_cmpd_name, y = 1, fill = Class.I)) +
  geom_tile() + 
  scale_fill_manual(values = class_colors) +
  labs(x = NULL, y = "") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14, color = "black"), 
    axis.line.y = element_blank(),
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = ggplot2::margin(t = 0, r = 5, b = 5, l = 5),
    legend.spacing.x = unit(2, "mm"),
    legend.key.size = unit(5, "mm"),
    legend.text = element_text(size = 14),
    legend.position = "bottom"
  )+
  guides(fill = guide_legend(
    title = NULL, 
    nrow = 1, 
    byrow = TRUE,
    label.position = "right" 
  ))


final_plot_1 <- p_bar / p_strip + plot_layout(heights = c(15, 1))



quartz(type = "pdf", file = "Anti_graph/fig4d_micro_class3_725_Academic_Standard_withoutwk_adj.pdf", width = 10, height = 6.5)


dev.off()












df_micro_pie <- df_annotation %>%
  group_by(Class.I) %>%
  summarise(count = n()) %>%
  mutate(prop = count / sum(count)) %>%
  arrange(desc(Class.I)) # 必须排序


df_micro_pie$ymax <- cumsum(df_micro_pie$prop)
df_micro_pie$ymin <- c(0, head(df_micro_pie$ymax, n=-1))

df_micro_pie$labelPosition <- (df_micro_pie$ymax + df_micro_pie$ymin) / 2

df_micro_pie$label <- paste0(df_micro_pie$Class.I, "\n", round(df_micro_pie$prop * 100, 1), "%")

df_micro_pie$type <- ifelse(df_micro_pie$prop < 0.04, "outside", "inside")

p1_fixed <- ggplot(df_micro_pie) +

  geom_rect(aes(ymin = ymin, ymax = ymax, xmin = 2.5, xmax = 4, fill = Class.I), 
            color = "white") +

  geom_text(
    data = subset(df_micro_pie, type == "inside"),
    aes(x = 3.2, y = labelPosition, label = label), # x=3.5 位于环中间
    size = 5, 
    color = "black", 
    fontface = "plain" # 【修改点】强制常规体 (非粗体)
  ) +
  
  geom_segment(
    data = subset(df_micro_pie, type == "outside"),
    aes(x = 4, xend = 4.2,             # 从环外沿(4)延伸到(4.2)
        y = labelPosition, yend = labelPosition), # 角度不变，径向延伸
    color = "black", 
    size = 0.5
  ) +
  
  geom_text(
    data = subset(df_micro_pie, type == "outside"),
    aes(x = 4.5, y = labelPosition, label = label), # x=4.25 在线头外面一点
    size = 5, 
    color = "black", 
    fontface = "plain", # 【修改点】强制常规体
    hjust = 0 # 稍微左对齐一点，防止撞线
  ) +
  
  scale_fill_manual(values = class_colors) +
  coord_polar(theta = "y") +

  xlim(c(1, 4.8)) + 
  
  theme_void() +
  theme(
    legend.position = "none", # 隐藏图例
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  ) 

print(p1_fixed)


ggsave('Anti_graph/fig4d_micro_circo.pdf', p1_fixed, width = 4, height = 4)







data_info = read.xlsx("../data/metabolites_info_624_THSBC_merge_corrected.xlsx")

data_medication = read.xlsx("Anti_results/Medication_Explaining_Metabolites_Ridge_withouwk_adj.xlsx")

data_clinical = read.xlsx("Anti_results/Clinical_Explaining_Metabolites_Ridge_withouwk_adj.xlsx")

data_lifestyle = read.xlsx("Anti_results/Lifestyle_Explaining_Metabolites_Ridge_withouwk_adj.xlsx")

data_baseline = read.xlsx("Anti_results/Baseline_Explaining_Metabolites_Ridge_withouwk_adj.xlsx")

data_dietary = read.xlsx("Anti_results/Dietary_Explaining_Metabolites_Ridge_withouwk_adj.xlsx")




info_clean <- data_info %>% 
  select(Index, sugg_cmpd_name, Class.I) %>% 
  distinct(Index, .keep_all = TRUE)

data_medication_annotated <- data_medication %>% 
  left_join(info_clean, by = c("metabolite" = "Index")) %>% 
  mutate(source = 'Medication use')

data_clinical_annotated <- data_clinical %>% 
  left_join(info_clean, by = c("metabolite" = "Index")) %>% 
  mutate(source = 'Clinical history')

data_lifestyle_annotated <- data_lifestyle %>% 
  left_join(info_clean, by = c("metabolite" = "Index")) %>% 
  mutate(source = 'Lifestyle factors')

data_baseline_annotated <- data_baseline %>% 
  left_join(info_clean, by = c("metabolite" = "Index")) %>% 
  mutate(source = 'Maternal baseline characteristics')

data_dietary_annotated <- data_dietary %>% 
  left_join(info_clean, by = c("metabolite" = "Index")) %>% 
  mutate(source = 'Dietary intake')

medication_all <- data_medication_annotated %>% 
  filter(period == 'All')

clinical_all <- data_clinical_annotated %>% 
  filter(period == 'All')

lifestyle_all <- data_lifestyle_annotated %>% 
  filter(period == 'All')

baseline_all <- data_baseline_annotated %>% 
  filter(period == 'All')

dietary_all <- data_dietary_annotated %>% 
  filter(period == 'All')


df_all <- rbind(medication_all, clinical_all, lifestyle_all, baseline_all, dietary_all)






df_all$R2 <- df_all$R2 * 100

df_all$sugg_cmpd_name <- str_replace_all(df_all$sugg_cmpd_name, regex(" Acid", ignore_case = TRUE), " acid")

name_mapping <- c(
  'Glycerophospho-N-Arachidonoyl Ethanolamine' = 'GP-NAE',
  'trans-resveratrol-3-O-sulfate' = 'Trans-resveratrol-3-O-sulfate',
  'Phosphatidylethanolamine lyso alkenyl 16:0' = 'LPE(P-16:0)',
  'Glycine deoxycholic acid'= 'GDCA',
  'Glycochenodeoxycholic acid'= 'GCDCA',
  'Deoxycholic acid'= 'DCA',
  'Taurocholic acid'= 'TCA',
  'Taurochenodesoxycholic acid'= 'TCDCA',
  'Glycoursodeoxycholic acid'= 'GUDCA',
  'Glycohyodeoxycholic acid'= 'GHDCA'
)

df_all$sugg_cmpd_name <- ifelse(
  df_all$sugg_cmpd_name %in% names(name_mapping),
  name_mapping[df_all$sugg_cmpd_name],
  df_all$sugg_cmpd_name
)



col_source <- c("Medication use" = "#687f99", "Clinical history" = "#9fc07f", 'Lifestyle factors' = '#e0756e', 
                'Maternal baseline characteristics' = '#e9ac70', 'Dietary intake' = '#9c91b8')

class_colors <- c("#FFFFB3", "#8DD3C7", "#BEBADA", "#FB8072", 
                  "#80B1D3", "#FDB462", "#FCCDE5", "#D9D9D9", 
                  "#BC80BD", "#CCEBC5", "#3498DB", "#FEB29B")
class_order <- c("AlAm", "AAD", "BA", "BSD", "CHD", "FA", 
                 "GP", "HRC", "NUC", "OAD", "SPL", "Others")

names(class_colors) <- class_order

stats_df <- df_all %>%
  group_by(sugg_cmpd_name) %>%
  summarise(total_R2 = sum(R2, na.rm = TRUE)) %>%
  # filter(total_R2 > 2) %>%
  arrange(desc(total_R2)) %>%
  slice_head(n = 30)
  


df_plot <- df_all %>%
  filter(sugg_cmpd_name %in% stats_df$sugg_cmpd_name)

df_plot$sugg_cmpd_name <- factor(df_plot$sugg_cmpd_name, 
                                 levels = stats_df$sugg_cmpd_name)

df_plot$source <- factor(df_plot$source, levels = c(
  "Maternal baseline characteristics",
  "Lifestyle factors",
  "Medication use",
  "Clinical history",
  "Dietary intake"
))


p_bar <- ggplot(df_plot, aes(x = sugg_cmpd_name, y = R2, fill = source)) +
  geom_col(width = 0.8, color = NA) + 
  scale_fill_manual(
    values = col_source,
    labels = c(
      "Maternal baseline characteristics" = "Sociodemographic and obstetrical characteristics", 
      "Lifestyle factors" = "Lifestyle factors",
      "Dietary intake" = "Dietary intake",
      "Clinical history" = "Clinical history",
      "Medication use" = "Medication use"
    )
  ) +
  # scale_fill_manual(values = col_source) +
  labs(y = "Variance explained (%)", x = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 4.7)) +
  theme_classic() +
  theme(
    # axis.title.y = element_text(size = 18, color = "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(), 
    axis.line.x = element_blank(),  
    legend.position = c(0.8, 0.85),
    legend.text = element_text(size = 11),
    # legend.title = element_blank(),
    axis.ticks.length.x = unit(0, "mm"), 
    plot.margin = margin(t = 5, r = 5, b = -5, l = 5)
  )


df_annotation <- df_plot %>%
  select(sugg_cmpd_name, Class.I) %>%
  distinct()

p_strip <- ggplot(df_annotation, aes(x = sugg_cmpd_name, y = 1, fill = Class.I)) +
  geom_tile() +
  # geom_tile(show.legend = FALSE) +
  scale_fill_manual(values = class_colors) +
  labs(x = NULL, y = "") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16, color = "black"), 
    axis.line.y = element_blank(),
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(t = -5, r = 5, b = 5, l = 5),
    legend.spacing.x = unit(2, "mm"),
    legend.key.size = unit(5, "mm"),
    legend.text = element_text(size = 14),
    legend.position = "bottom"
  )+
  guides(fill = guide_legend(
    title = NULL,
    nrow = 1,
    byrow = TRUE,
    label.position = "right"
  ))


final_plot_2 <- p_bar / p_strip + plot_layout(heights = c(15, 1))


print(final_plot_2)




ggsave('Anti_graph/fig4d_factors_class3_withouwk_adj.pdf', final_plot_2, width = 10, height = 6.2)




combined_plot <- final_plot_1 | final_plot_2 + 
  plot_annotation(tag_levels = 'A') 

# 打印查看效果
print(combined_plot)

# 保存最终的拼接大图（注意高度 height 翻倍了）
quartz(type = "pdf", file = "Anti_graph/fig4_combined_factors_micro_adj.pdf", width = 20, height = 6.7)
print(combined_plot)
dev.off()




df_micro_pie <- df_annotation %>%
  group_by(Class.I) %>%
  summarise(count = n()) %>%
  mutate(prop = count / sum(count)) %>%
  arrange(desc(Class.I)) # 必须排序


df_micro_pie$ymax <- cumsum(df_micro_pie$prop)
df_micro_pie$ymin <- c(0, head(df_micro_pie$ymax, n=-1))

df_micro_pie$labelPosition <- (df_micro_pie$ymax + df_micro_pie$ymin) / 2

df_micro_pie$label <- paste0(df_micro_pie$Class.I, "\n", round(df_micro_pie$prop * 100, 1), "%")

df_micro_pie$type <- ifelse(df_micro_pie$prop < 0.04, "outside", "inside")

p1_fixed <- ggplot(df_micro_pie) +
  
  geom_rect(aes(ymin = ymin, ymax = ymax, xmin = 2.5, xmax = 4, fill = Class.I), 
            color = "white") +
  
  geom_text(
    data = subset(df_micro_pie, type == "inside"),
    aes(x = 3.2, y = labelPosition, label = label), # x=3.5 位于环中间
    size = 5, 
    color = "black", 
    fontface = "plain" # 【修改点】强制常规体 (非粗体)
  ) +
  
  geom_segment(
    data = subset(df_micro_pie, type == "outside"),
    aes(x = 4, xend = 4.2,             # 从环外沿(4)延伸到(4.2)
        y = labelPosition, yend = labelPosition), # 角度不变，径向延伸
    color = "black", 
    size = 0.5
  ) +
  
  geom_text(
    data = subset(df_micro_pie, type == "outside"),
    aes(x = 4.5, y = labelPosition, label = label), # x=4.25 在线头外面一点
    size = 5, 
    color = "black", 
    fontface = "plain", # 【修改点】强制常规体
    hjust = 0 # 稍微左对齐一点，防止撞线
  ) +
  
  scale_fill_manual(values = class_colors) +
  coord_polar(theta = "y") +
  
  xlim(c(1, 4.8)) + 
  
  theme_void() +
  theme(
    legend.position = "none", # 隐藏图例
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  ) 

print(p1_fixed)


ggsave('Anti_graph/fig4d_factors_circo.pdf', p1_fixed, width = 4, height = 4)




