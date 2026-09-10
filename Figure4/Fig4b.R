rm(list = ls())
setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code')

library(dplyr)
library(stringr)
library(reshape2)
library(openxlsx)
library(ggplot2)
library(survey)
library(lme4)
library(lmerTest)
library(patchwork)
library(UpSetR)
library(tidyr)

df_t1 = read.xlsx('../Anti_data/preprocessed/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov_class3.xlsx', sheet = 7) %>% mutate(period = 'T1')
df_t2 = read.xlsx('../Anti_data/preprocessed/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov_class3.xlsx', sheet = 8) %>% mutate(period = 'T2')
df_t3 = read.xlsx('../Anti_data/preprocessed/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov_class3.xlsx', sheet = 9) %>% mutate(period = 'T3')
df_pool = read.xlsx('../Anti_data/preprocessed/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov_class3.xlsx', sheet = 14) %>% mutate(period = 'Pooled')

df_t1$inf_flag <- NULL
df_t2$inf_flag <- NULL
df_t3$inf_flag <- NULL
df = rbind(df_t1, df_t2, df_t3, df_pool)

significant_compounds <- df %>%
  filter(pv_adj_FDR < 0.05) %>%
  group_by(sugg_cmpd_name) %>%
  summarise(n_periods = n_distinct(period)) %>%
  filter(n_periods == 4) %>%
  pull(sugg_cmpd_name)

significant_compounds

cut.pv = 0.05
cut.beta = 0.0 
df = df %>% mutate(col = case_when(pv_adj_FDR < cut.pv & beta > cut.beta ~ 'Elevated',
                                   pv_adj_FDR < cut.pv & beta < -cut.beta ~ 'Reduced',
                                   .default = 'notsig'))

df_pos = df %>% filter(col=='Elevated')
df_pos = split(df_pos$meta_name, df_pos$period)
pdf("Anti_graph/upset_meta_up_GH_class3.pdf", width = 7.05, height = 6.49)
bar_col = "darkred"
upset(fromList(df_pos), sets = names(df_pos), 
      keep.order = T, decreasing = TRUE,
      set_size.show = T, set_size.scale_max = 60,
      matrix.color = bar_col, sets.bar.color = bar_col,
      main.bar.color = bar_col, shade.color = bar_col,
      shade.alpha = 0.15, text.scale = c(3, 3, 3, 3, 3, 4),
      point.size = 4, line.size = 0.9,
      order.by = "freq", mb.ratio = c(0.7,0.3), 
      mainbar.y.label = "Intersection size", sets.x.label = "Set size")

while (!is.null(dev.list())) dev.off()

df_neg = df %>% filter(col=='Reduced')
df_neg = split(df_neg$meta_name, df_neg$period)
pdf("Anti_graph/upset_meta_down_GH_class3.pdf", width = 7.05, height = 6.49)
bar_col = "#4682B4"
upset(fromList(df_neg), sets = names(df_neg), 
      keep.order = T, decreasing = TRUE,
      set_size.show = T, set_size.scale_max = 15,
      matrix.color = bar_col, sets.bar.color = bar_col,
      main.bar.color = bar_col, shade.color = bar_col,
      shade.alpha = 0.15, text.scale = c(3, 3, 3, 3, 3, 4),
      point.size = 3.5, line.size = 0.9,
      order.by = "freq", mb.ratio = c(0.7,0.3), 
      mainbar.y.label = "Intersection size", sets.x.label = "Set size")

while (!is.null(dev.list())) dev.off()

df_diff = df %>% filter(col!='notsig') %>% 
  select(period, Class.I) %>% 
  group_by(period, Class.I) %>% 
  summarise(count = n(), .groups = 'drop')

class_o = df_diff %>% group_by(Class.I) %>% 
  summarise(mean_count = mean(count), .groups = 'drop') %>% 
  arrange(desc(mean_count))
other_class = class_o %>% filter(mean_count<2) %>% pull(Class.I)

df_pos = df %>% filter(col=='Elevated') %>% 
  mutate(Class.I = ifelse(Class.I %in% other_class, 'Others', Class.I)) %>% 
  group_by(period, Class.I) %>% 
  summarise(count = n(), .groups = 'drop')

df_neg = df %>% filter(col=='Reduced') %>% 
  mutate(Class.I = ifelse(Class.I %in% other_class, 'Others', Class.I)) %>% 
  group_by(period, Class.I) %>% 
  summarise(count = n(), .groups = 'drop')


colorset = c( "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#3498DB", "#FEB29B")
order = c( "AAD", "BA", "BSD", "CHD", "FA", "GP", "HRC", "NUC", "OAD", "SPL", "Others")

color_mapping = setNames(colorset, order)

#########

all_classes = unique(c(df_pos$Class.I, df_neg$Class.I)) # Get all categories
all_classes = c(all_classes, "Others") # Ensure the "Others" category is included
all_classes = unique(all_classes)

df_pos = df_pos %>%
  ungroup() %>% 
  complete(period, Class.I = all_classes, fill = list(count = 0)) 

df_neg = df_neg %>%
  ungroup() %>% 
  complete(period, Class.I = all_classes, fill = list(count = 0))

ordered_classes = c(
  order[order %in% setdiff(all_classes, c("BA", "Others"))],
  "BA", 
  "Others"
)

df_pos$Class.I = factor(df_pos$Class.I, levels = ordered_classes)
df_neg$Class.I = factor(df_neg$Class.I, levels = ordered_classes)

class_col = color_mapping[ordered_classes]

p_pos = ggplot(df_pos, aes(x=period, y = count, fill = Class.I))+
  geom_bar(stat="identity", position='stack', width = 0.7)+
  scale_fill_manual(limits = all_classes, values = class_col)+
  labs(x='', y='The number of metabolites', fill='Class', title = 'Elevated')+
  theme_bw() +theme(panel.grid= element_blank(),
                    axis.text.y = element_text(size = 18, color = 'black'),
                    axis.text.x = element_text(size = 16, angle = 35, color = 'black', margin = ggplot2::margin(t = 8)),
                    axis.title.y = element_text(size = 18, color = 'black'))

p_neg = ggplot(df_neg, aes(x=period, y = count, fill = Class.I))+
  geom_bar(stat="identity", position='stack', width = 0.7)+
  scale_fill_manual(limits = all_classes, values = class_col)+
  labs(x='', y='', fill='Class', title = 'Reduced')+
  theme_bw() +theme(panel.grid= element_blank(),
                    axis.text.y = element_text(size = 18, color = 'black'),
                    axis.text.x = element_text(size = 16, angle = 35, color = 'black', margin = ggplot2::margin(t = 8)))

combined_plot = (p_pos / p_neg) + 
  plot_layout(guides = 'collect') & 
  theme(
    legend.text = element_text(size = 16), 
    legend.title = element_text(size = 16),
    plot.margin = ggplot2::margin(5, 5, 5, 5),
    plot.title = element_text(size = 16, hjust = 0.5)
  )

ggsave('Anti_graph/diff_metabolite_class_logz_GH_class3_725.pdf', combined_plot, width = 4.8, height = 7.)