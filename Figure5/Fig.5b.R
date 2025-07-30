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



df_t1 = read.xlsx('../data/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov.xlsx', sheet = 7) %>% mutate(period = 'T1')
df_t2 = read.xlsx('../data/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov.xlsx', sheet = 8) %>% mutate(period = 'T2')
df_t3 = read.xlsx('../data/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov.xlsx', sheet = 9) %>% mutate(period = 'T3')
# sheet 1/2/3 for HDP 4/5/6 for PE 7/8/9 for GH

df = rbind(df_t1, df_t2, df_t3)
significant_compounds <- df %>%
  filter(pv_adj_FDR < 0.05) %>%
  group_by(sugg_cmpd_name) %>%
  summarise(periods = list(unique(period)), count = n()) %>%
  filter(all(c("T1", "T2", "T3") %in% unlist(periods)) & count >= 3) %>%
  pull(sugg_cmpd_name)


significant_compounds


cut.pv = 0.05
cut.beta = 0.0 
df = df %>% mutate(col = case_when(pv_adj_FDR < cut.pv & beta > cut.beta ~ 'up',
                                   pv_adj_FDR < cut.pv & beta < -cut.beta ~ 'down',
                                   .default = 'notsig'))


library(UpSetR)
df_pos = df %>% filter(col=='up')
df_pos = split(df_pos$meta_name, df_pos$period)
pdf("graph/upset_meta_up_GH.pdf", width = 7.05, height = 5.49)
bar_col = "darkred"
upset(fromList(df_pos), sets = names(df_pos), 
      keep.order = T, decreasing = TRUE,
      set_size.show = T, set_size.scale_max = 150,
      matrix.color = bar_col, sets.bar.color = bar_col,
      main.bar.color = bar_col, shade.color = bar_col,
      shade.alpha = 0.15, text.scale = c(3, 3, 3, 3, 3, 4),
      point.size = 4, line.size = 0.9,
      order.by = "freq", mb.ratio = c(0.7,0.3), 
      mainbar.y.label = "Intersection size", sets.x.label = "Set size")

while (!is.null(dev.list())) dev.off()

df_neg = df %>% filter(col=='down')
df_neg = split(df_neg$meta_name, df_neg$period)
pdf("Anti_graph/upset_meta_down_GH.pdf", width = 6.05, height = 5.49)
bar_col = "#4682B4"
upset(fromList(df_neg), sets = names(df_neg), 
      keep.order = T, decreasing = TRUE,
      set_size.show = T, set_size.scale_max = 20,
      matrix.color = bar_col, sets.bar.color = bar_col,
      main.bar.color = bar_col, shade.color = bar_col,
      shade.alpha = 0.15, text.scale = c(3, 3, 3, 3, 3, 4),
      point.size = 3.5, line.size = 0.9,
      order.by = "freq", mb.ratio = c(0.7,0.3), 
      mainbar.y.label = "Intersection size", sets.x.label = "Set size")


while (!is.null(dev.list())) dev.off()


library(tidyr)
df_diff = df %>% filter(col!='notsig') %>% 
  select(period, Class.I) %>% 
  group_by(period, Class.I) %>% 
  summarise(count = n())
class_o = df_diff %>% group_by(Class.I) %>% 
  summarise(mean_count = mean(count)) %>% 
  arrange(desc(mean_count))
other_class = class_o %>% filter(mean_count<2) %>% pull(Class.I)

df_pos = df %>% filter(col=='up') %>% 
  mutate(Class.I = ifelse(Class.I %in% other_class, 'Others', Class.I)) %>% 
  group_by(period, Class.I) %>% 
  summarise(count = n())
df_neg = df %>% filter(col=='down') %>% 
  mutate(Class.I = ifelse(Class.I %in% other_class, 'Others', Class.I)) %>% 
  group_by(period, Class.I) %>% 
  summarise(count = n())

colorset = c("#FFFFB3", "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#3498DB", "#FEB29B")#, "#FFD19D"
order = c("AA", "AAM", "BA", "BSD", "CHM", "FA", "GP", "HetC", "HHC", "NTM", "OAD", "SPL", "Others")

color_mapping = setNames(colorset, order)


#########

all_classes = unique(c(df_pos$Class.I, df_neg$Class.I))
all_classes = c(all_classes, "Others") 
all_classes = unique(all_classes)

df_pos = df_pos %>%
  ungroup() %>%
  complete(period, Class.I = all_classes, fill = list(count = 0))


df_neg = df_neg %>%
  ungroup() %>% 
  complete(period, Class.I = all_classes, fill = list(count = 0)) 


all_classes = unique(c(df_pos$Class.I, df_neg$Class.I))


ordered_classes = c(
  order[order %in% setdiff(all_classes, c("BA", "Others"))],
  "BA", 
  "Others"
)

df_pos$Class.I = factor(df_pos$Class.I, levels = ordered_classes)
df_neg$Class.I = factor(df_neg$Class.I, levels = ordered_classes)


class_col = color_mapping[ordered_classes]



#########

p_pos = ggplot(df_pos, aes(x=period, y = count, fill = Class.I))+
  geom_bar(stat="identity", position='stack', width = 0.7)+
  scale_fill_manual(limits = all_classes, values = class_col)+
  labs(x='', y='The number of metabolites', fill='Class', title = 'Up')+
  theme_bw() +theme(panel.grid= element_blank(),
                    axis.text.y = element_text(size = 14, color = 'black'),
                    axis.text.x = element_text(size = 14, color = 'black'))

p_neg = ggplot(df_neg, aes(x=period, y = count, fill = Class.I))+
  geom_bar(stat="identity", position='stack', width = 0.7)+
  scale_fill_manual(limits = all_classes, values = class_col)+
  labs(x='', y='', fill='Class', title = 'Down')+
  theme_bw() +theme(panel.grid= element_blank(),
                    axis.text.y = element_text(size = 14, color = 'black'),
                    axis.text.x = element_text(size = 14, color = 'black'))


combined_plot = p_pos + p_neg + plot_layout(guides = 'collect')



