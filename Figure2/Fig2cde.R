#######new data ############


rm(list = ls())
library(tidyverse)
library(data.table)
library(readxl)
library(writexl)
library(haven)

library(ggplot2)
library(patchwork)
library(ggsignif)
library(gghalves)
library(ggh4x)
library(vegan)
library(compositions)

library(lme4)
library(lmerTest)
library(cowplot)
library(purrr)

setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code')
source('../R_code/myUtils.R')
set.seed(123)

bcdist = read.csv('data/Weighted_ALL_ITS.csv')
# bcdist = read.csv('data/Weighted_ALL_Taxonomy.csv')
bcdist <- bcdist[!is.na(bcdist$wt1), ]


bcdist$period <- sub("^V", "T", as.character(bcdist$period))
bcdist <- bcdist %>% 
  # 1. 删除 HDP 为 NA 的行
  filter(!is.na(HDP)) %>% 
  # 2. 按规则生成 group 列
  mutate(
    group = case_when(
      HDP == 1 & preeclampsia == 1 ~ "PE",
      HDP == 1 & preeclampsia == 0 ~ "GH",
      HDP == 0                    ~ "NP",
      TRUE                        ~ NA_character_   # 兜底情况（以防出现其他取值）
    )
  )

# table(bcdist$group)

load('data/beta_diversity_rawdata_forplot.RData')

# species_col  <- c(grep("^s\\d+$", names(bcdist), value = TRUE))
# abund_tab = decostand(select(bcdist, all_of(species_col)), 'total', 1)

# abund_tab = decostand(select(bcdist, all_of(bact_mgs_s_uni)), 'total', 1)

abund_tab = decostand(select(bcdist, all_of(fung_its_uni)), 'total', 1)

eucdist <- vegdist(abund_tab, method = "bray")



# load('data/weighted-unifrac-distance-matrix.rdata')





# keep_cols <- c("id", "period", "group", all_of(bact_mgs_s_uni))
keep_cols <- c("id", "period", "group", all_of(fung_its_uni))

species_sample <- bcdist[, keep_cols, drop = FALSE]

# 假设检验
#MGS
# Adonis R2 for period = 0.054, P = 0.001
# Adonis R2 for group = 0.004, P = 0.001

#MGG
# Adonis R2 for period = 0.073, P = 0.001
# Adonis R2 for group = 0.005, P = 0.001

# ITS
# "Adonis R2 for period = 0.004, P = 0.001"
# "Adonis R2 for Group = 0.002, P = 0.004"
period_test = adonis2(eucdist~species_sample$period, permutations = 999, parallel = 8)
R2_time = period_test$R2[1]; pv_time = period_test$`Pr(>F)`[1]
group_test = adonis2(eucdist~species_sample$group, permutations = 999, parallel = 8)
R2_group = group_test$R2[1]; pv_group = group_test$`Pr(>F)`[1]
sprintf('Adonis R2 for period = %.3f, P = %.3f', R2_time, pv_time)
sprintf('Adonis R2 for Group = %.3f, P = %.3f', R2_group, pv_group)

# 时期 分别检验
# Adonis R2 for ICP = 0.000, P = 0.025
# Adonis R2 for ICP = 0.004, P = 0.001
# Adonis R2 for ICP = 0.006, P = 0.001

# MGS
# Adonis R2 for Group = 0.011, P = 0.001
# Adonis R2 for Group = 0.003, P = 0.792
# Adonis R2 for Group = 0.004, P = 0.247

# MGG
# Adonis R2 for Group = 0.014, P = 0.001
# Adonis R2 for Group = 0.002, P = 0.891
# Adonis R2 for Group = 0.005, P = 0.189

# ITS
# Adonis R2 for Group = 0.005, P = 0.001
# Adonis R2 for Group = 0.002, P = 0.233
# Adonis R2 for Group = 0.003, P = 0.170
# species_sample$period <- sub("^V", "T", as.character(species_sample$period))
for (t in c('T1','T2','T3')){
  set.seed(1111)
  idx_t = which(species_sample$period == t)
  dist_t = as.dist(as.matrix(eucdist)[idx_t, idx_t])
  group_t = species_sample$group[idx_t]
  group_test = adonis2(dist_t~group_t, permutations = 999, parallel = 8)
  R2_group = group_test$R2[1]; pv_group = group_test$`Pr(>F)`[1]
  cat(sprintf('Adonis R2 for Group = %.3f, P = %.3f\n', R2_group, pv_group))
}

# 疾病组 分别检验
# MGS
# Adonis R2 for period = 0.023, P = 0.038
# Adonis R2 for period = 0.031, P = 0.002
# Adonis R2 for period = 0.061, P = 0.001

# MGG
# Adonis R2 for period = 0.030, P = 0.008
# Adonis R2 for period = 0.035, P = 0.002
# Adonis R2 for period = 0.082, P = 0.001

# ITS
# Adonis R2 for period = 0.014, P = 0.412
# Adonis R2 for period = 0.018, P = 0.310
# Adonis R2 for period = 0.005, P = 0.001
for (case in c('PE',"GH", 'NP')){
  set.seed(1111)
  idx_t = which(species_sample$group == case)
  dist_t = as.dist(as.matrix(eucdist)[idx_t, idx_t])
  group_t = species_sample$period[idx_t]
  group_test = adonis2(dist_t~group_t, permutations = 999, parallel = 2)
  R2_group = group_test$R2[1]; pv_group = group_test$`Pr(>F)`[1]
  cat(sprintf('Adonis R2 for period = %.3f, P = %.3f\n', R2_group, pv_group))
}

# ptColors3 = c('#7a94ba','#f8c672','#B5665D')

pcoa_w_unifrac = cmdscale(eucdist, k = 2, eig = TRUE)
# save(pcoa_w_unifrac, file = "data/MGS_S_bray_pcoa.RData")

contr = pcoa_w_unifrac$eig/sum(pcoa_w_unifrac$eig) * 100
head(contr)

# load("data/MGS_S_bray_pcoa.RData")



# 绘图
col_group = c('PE'='#d73027',"GH" = "#f8c672",'NP'='#7a94ba')
col_period = shades::saturation(c('T1'='#f0d9a2','T2'='#a8d5c7','T3'='#c3dff5'), 0.6)
colnames(pcoa_w_unifrac$points) = paste0('PCOA',1:2)
df_pcoa = cbind(species_sample %>% select(id, period, group), pcoa_w_unifrac$points[,1:2])
# 主图
p_main = ggplot(df_pcoa, aes(x = PCOA1, y = PCOA2)) +
  geom_point(aes(color = period, shape = group), size = 1.5, alpha = 0.8)+
  scale_x_continuous(position = "top")+
  annotate('text', x=-0.5, y=0.5, size=5, hjust = 0,
           label=expression(paste(R^2,'= 0.002, ',italic('P'),'= 0.004 for group')))+
  annotate('text', x=-0.5, y=0.45, size=5, hjust = 0,
           label=expression(paste(R^2,'= 0.004, ',italic('P'),'< 0.001 for period')))+
  scale_color_manual(values = col_period, name = "Study", guide = "none")+
  labs(x=sprintf('PCOA1 (%.1f%%)', contr[1]),
       y=sprintf('PCOA2 (%.1f%%)', contr[2])) +
  theme_test(base_size = 15) + 
  theme(
    legend.position = c(0.15, 0.1),
    legend.background = element_blank(), legend.title = element_blank(),
    axis.text = element_blank(), axis.ticks = element_blank(),
    plot.margin = margin(3,3,3,3)
  ) 
# 右侧1
df_pcoa$group <- factor(df_pcoa$group, levels = c("NP", "GH", "PE"))
df_pcoa$period <- factor(df_pcoa$period, levels = c("T1", "T2", "T3"))
r1 = ggplot(df_pcoa, aes(x = group, y = PCOA2, fill = group, color = group)) +
  geom_boxplot(size= 0.8, width = 0.5, alpha = 0.1) +
  scale_fill_manual(values = col_group)+
  scale_color_manual(values = col_group)+
  theme_test(base_size = 15)+
  theme(axis.title = element_blank(),
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.margin = margin(0,3.2,0,0))
# 右侧2
r2 = ggplot(df_pcoa, aes(x = period, y = PCOA2, fill = period, color = period)) +
  geom_boxplot(size= 0.8, width = 0.5, alpha = 0.1)+
  scale_y_continuous(position = "right")+
  scale_fill_manual(values = col_period)+
  scale_color_manual(values = col_period)+
  theme_test(base_size = 15)+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(0,0,0,0))

p_right = (r1 | r2) + plot_layout(widths = c(3, 3))
# 下侧1
b1 = ggplot(df_pcoa, aes(x = group, y = PCOA1, fill = group, color = group)) +
  geom_boxplot(size= 0.8, width = 0.5, alpha = 0.1)+
  scale_x_discrete(limits = c('NP', 'GH', "PE"))+
  scale_fill_manual(values = col_group)+
  scale_color_manual(values = col_group)+
  coord_flip()+
  theme_test(base_size = 15)+
  labs(fill = 'Group', color = 'Group')+
  theme(axis.title = element_blank(),
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.box.just = 'top', legend.justification = c(0,1),
        plot.margin = margin(0,0,0,0))
# 下测2
b2 = ggplot(df_pcoa, aes(x = period, y = PCOA1, fill = period, color = period)) +
  geom_boxplot(size= 0.8, width = 0.5, alpha = 0.1) +
  scale_x_discrete(limits = c('T3', 'T2', 'T1'))+
  scale_fill_manual(values = col_period)+
  scale_color_manual(values = col_period)+
  coord_flip()+
  theme_test(base_size = 15)+
  labs(fill = 'Period', color = 'Period')+
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        legend.box.just = 'top', legend.justification = c(0,1),
        plot.margin = margin(0,0,0,0))

p_bottom = (b1 / b2) & theme(legend.position = 'none')
legend_bottom <- cowplot::plot_grid(get_legend(b1), get_legend(b2), nrow = 1, align = 'v', rel_widths = c(1, 0.6, 1))

(p_main + p_right + p_bottom + legend_bottom) +
  plot_layout(heights = c(2.7, 0.8), widths = c(2.7,1))
ggsave('Anti_graph/fig2_bray_PCOA_ITS.pdf', width = 6.2, height = 6.2)









library(rstatix) # Great for statistical tables

# 1. Check Global Significance (Is there any difference among the 3 groups?)
# For PCoA1
res_kruskal_pcoa1 <- df_pcoa %>% kruskal_test(PCOA1 ~ group)
print(res_kruskal_pcoa1)

# For PCoA2
res_kruskal_pcoa2 <- df_pcoa %>% kruskal_test(PCOA2 ~ group)
print(res_kruskal_pcoa2)

# 2. Pairwise Comparisons (Which specific groups are different?)
# This will give you p-values for NP vs GH, GH vs PE, etc.
stat_pcoa1 <- df_pcoa %>% 
  wilcox_test(PCOA1 ~ group) %>%
  adjust_pvalue(method = "BH") %>% # Benjamini-Hochberg correction
  add_significance()

stat_pcoa2 <- df_pcoa %>% 
  wilcox_test(PCOA2 ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# View results
print(stat_pcoa1)
print(stat_pcoa2)

library(ggpubr)




















rm(list = ls())
library(tidyverse)
library(data.table)
library(readxl)
library(writexl)
library(haven)

library(ggplot2)
library(patchwork)
library(ggsignif)
library(gghalves)
library(ggh4x)
library(vegan)
library(compositions)

library(lme4)
library(lmerTest)
library(cowplot)
library(purrr)
library(rstatix)   # 显著性检验用
library(boot)      # Bootstrap用
library(ggpubr)    # 用于 stat_pvalue_manual 和 stat_compare_means

setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code')
source('../R_code/myUtils.R')
source('help_func.R') # 提前加载 help_func
set.seed(123)

# ==============================================================================
# 1. 数据读取与预处理
# ==============================================================================
# bcdist = read.csv('data/Weighted_ALL_ITS.csv')
bcdist = read.csv('data/Weighted_ALL_Taxonomy.csv')
bcdist <- bcdist[!is.na(bcdist$wt2), ]
bcdist$period <- sub("^V", "T", as.character(bcdist$period))

bcdist <- bcdist %>% 
  filter(!is.na(HDP)) %>% 
  mutate(
    group = case_when(
      HDP == 1 & preeclampsia == 1 ~ "PE",
      HDP == 1 & preeclampsia == 0 ~ "GH",
      HDP == 0                     ~ "NP",
      TRUE                         ~ NA_character_   
    )
  )

load('data/beta_diversity_rawdata_forplot.RData')

abund_tab = decostand(select(bcdist, all_of(bact_mgs_s_uni)), 'total', 1)
# abund_tab = decostand(select(bcdist, all_of(fung_its_uni)), 'total', 1)
eucdist <- vegdist(abund_tab, method = "bray")

keep_cols <- c("id", "period", "group", all_of(bact_mgs_s_uni))
species_sample <- bcdist[, keep_cols, drop = FALSE]
species_sample$group <- factor(species_sample$group, levels = c("NP", "GH", "PE"))


load("data/MGS_S_bray_pcoa.RData")
# load("data/ITS_bray_pcoa.RData")
pcoa1_per <- round(pcoa_w_unifrac$eig[1] / sum(pcoa_w_unifrac$eig) * 100, 2) 
ylab_text <- paste("PCoA1 (", pcoa1_per, "%)", sep = "")

colnames(pcoa_w_unifrac$points) = paste0('PCOA',1:2)
df_pcoa = cbind(species_sample %>% select(id, period, group), pcoa_w_unifrac$points[,1:2])


# ==============================================================================
# 2. 绘制第一张图 (p_pcoa1_box) - 按时间点分组
# ==============================================================================
df_pcoa$group <- factor(df_pcoa$group, levels = c("NP", "GH", "PE"))
col_group = colorspace::lighten(c('#7a94ba',"#f8c672",'#d73027'), amount = 0.1)

# 使用 rstatix 进行各时间点内的两两 Wilcoxon 检验，并进行 BH 校正
stat_test_pcoa <- df_pcoa %>%
  group_by(period) %>%
  wilcox_test(PCOA1 ~ group, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "period", dodge = 0.8) # 计算箱线图 dodge 后的 P 值标记位置

p_pcoa1_box <- ggplot(df_pcoa, aes(x = period, y = PCOA1)) +
  geom_boxplot(aes(fill = group, color = group), 
               size = 0.8, width = 0.6, alpha = 0.1, 
               position = position_dodge(0.8)) +
  scale_fill_manual(values = col_group) +
  scale_color_manual(values = col_group) +
  # 添加显著性标记（隐藏不显著的比较线让图形更干净）
  stat_pvalue_manual(stat_test_pcoa, label = "p.adj.signif", size = 4.5,
                     tip.length = 0.01, step.increase = 0.07, 
                     hide.ns = TRUE, bracket.size = 0.5) +
  labs(y = ylab_text, x = "Period") + 
  theme_classic(base_size = 12) +
  theme(
    axis.title.y = element_text(color = "black", size = 12), 
    axis.text.x  = element_text(color = "black", size = 12),
    axis.text.y  = element_text(color = "black", size = 12), 
    axis.title.x = element_blank(),
    legend.position = 'none',      # 保留第一张图的图例在顶部
    legend.title = element_blank()
  )


# ==============================================================================
# 3. 绘制第二张图 (p_trend) - 放在中间
# ==============================================================================
df_trend = data.frame() 
for (t in c('T1','T2','T3')){
  for (g in c('PE', 'GH', 'NP')){
    vdata_tmp = species_sample %>% filter(group==g & period==t)
    tem_abund = decostand(select(vdata_tmp, all_of(bact_mgs_s_uni)), 'total', 1)
    vdist = vegdist(tem_abund, 'bray')
    vdist = as.matrix(vdist)
    df_trend = rbind(df_trend, data.frame(period=t, group=g, dist=vdist[lower.tri(vdist)]))
  }
}

df_trend <- df_trend %>% mutate(group = factor(group, levels = c("NP", "GH", "PE")))

stat_test <- df_trend %>%
  group_by(period) %>%
  wilcox_test(dist ~ group, p.adjust.method = "BH") %>%
  add_significance("p.adj") 

sig_data <- stat_test %>%
  filter(p.adj.signif != "ns") %>% 
  mutate(
    y_pos = case_when(
      # group1 == "NP" & group2 == "PE" ~ 0.89,
      # group1 == "GH" & group2 == "PE" ~ 0.87,
      # group1 == "NP" & group2 == "GH" ~ 0.85,
      # TRUE ~ 0.85
      group1 == "NP" & group2 == "PE" ~ 0.78,
      group1 == "GH" & group2 == "PE" ~ 0.765,
      group1 == "NP" & group2 == "GH" ~ 0.75,
      TRUE ~ 0.75
    ),
    comp_color = case_when(
      group1 == "NP" & group2 == "PE" ~ "#c4827f", 
      group1 == "GH" & group2 == "PE" ~ "#cc79a7", 
      group1 == "NP" & group2 == "GH" ~ "#e8bf64", 
      TRUE ~ "black"
    )
  )

set.seed(1111)
df_med = df_trend %>% group_by(period, group) %>% 
  do(med_boot = boot.ci(boot(.$dist, median_func, R = 100), type='perc')) %>% 
  mutate(med = med_boot$t0, med_lci = med_boot$percent[4], med_uci = med_boot$percent[5]) %>%
  mutate(group = factor(group, levels = c("NP", "GH", "PE")))

p_trend <- ggplot(df_med, aes(x=period, fill=group, group=group))+
  geom_ribbon(aes(ymin=med_lci, ymax=med_uci), alpha=0.15)+
  geom_line(aes(y=med, color=group))+
  geom_point(aes(y=med, color=group, shape=group), size=1.2)+
  scale_x_discrete(expand = c(0,0.2))+
  # scale_y_continuous(limits = c(NA, 0.90))+ 
  scale_y_continuous(limits = c(NA, 0.78))+
  scale_color_manual(values = col_group)+
  scale_fill_manual(values = col_group)+
  scale_shape_manual(values = c("NP" = 16, "GH" = 15, "PE" = 17))+
  geom_text(data = sig_data,
            aes(x = period, y = y_pos, label = p.adj.signif, color = I(comp_color)),
            size = 4.5, inherit.aes = FALSE, fontface = "bold", vjust = 0.5) +
  labs(y = 'Bray-Curtis dissimilarity')+
  theme_classic(base_line_size = 0.5)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black", 
                                    vjust = 1.9, hjust = 0.5, angle = 90),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "grey10", 
                                   vjust = 0.5, hjust = 0.5, angle = 0),
        axis.ticks.length = unit(1.5, 'mm'),
        legend.position = 'none') # 去除图例


# ==============================================================================
# 4. 绘制第三张图 (p1) - 放在最右侧
# ==============================================================================
id_PE = unique(species_sample %>% filter(group=='PE') %>% pull(id))
id_GH = unique(species_sample %>% filter(group=='GH') %>% pull(id))
id_nHDP = unique(species_sample %>% filter(group=='NP') %>% pull(id))

vdist_matrix = as.matrix(eucdist)
within_id = get_dist_pair(species_sample$id)
dist_id = vdist_matrix[as.matrix(within_id[,1:2])]
dist_PE = vdist_matrix[as.matrix(within_id[within_id$group_row %in% id_PE, 1:2])]
dist_GH = vdist_matrix[as.matrix(within_id[within_id$group_row %in% id_GH, 1:2])]
dist_nHDP = vdist_matrix[as.matrix(within_id[within_id$group_row %in% id_nHDP, 1:2])]

within_period = get_dist_pair(species_sample$period)
dist_period = vdist_matrix[as.matrix(within_period[,1:2])]

df_stability = rbind(data.frame(type='Within\nindividual', dist=dist_id),
                     data.frame(type='Within PE\nindividual', dist=dist_PE),
                     data.frame(type='Within GH\nindividual', dist=dist_GH),
                     data.frame(type='Within NP\nindividual', dist=dist_nHDP),
                     data.frame(type='Within\nperiod', dist=dist_period))

p1 = ggplot(df_stability %>% filter(type %in% c('Within\nindividual','Within\nperiod')), 
            aes(x=type, y=dist, color=type))+
  geom_boxplot(size = 0.8, width = 0.5, alpha = 0.1, outliers = FALSE)+
  scale_y_continuous(expand = c(0, 0.02), limits = c(0, 1.08))+
  scale_color_manual(values = c('#1d2d44','#8daa9d'))+
  geom_signif(comparisons = list(c("Within\nindividual", "Within\nperiod")), 
              annotations = c('****'), y_position = 0.98,
              color = 'black', tip_length = 0.03, textsize = 4.5)+
  labs(y = 'Bray-Curtis dissimilarity')+
  theme_classic(base_line_size = 0.5)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black", 
                                    vjust = 1.9, hjust = 0.5, angle = 90),
        axis.text.x = element_text(size = 12, color = "black", angle = 0, vjust=0.65, hjust=0.5), 
        axis.text.y = element_text(size = 12, color = "grey10",
                                   vjust = 0.5, hjust = 0.5, angle = 0),
        axis.ticks.length = unit(1.5, 'mm'), 
        legend.position = 'none')


# ==============================================================================
# 5. 拼图及保存 (Patchwork)
# ==============================================================================
# 使用 plot_layout 调整三张图的宽度比例
p_plot = p_pcoa1_box + p_trend + p1 + plot_layout(widths = c(2.3, 2.3, 1.4))

# 调整了画板的尺寸以适应新的比例
ggsave('Anti_graph/fig2d_betadist3_MGS_horiz_merged_multi.pdf', p_plot, width = 8.5, height = 3.1)









###########################PCOA ITS MGS ###############

rm(list = ls())
library(tidyverse)
library(data.table)
library(readxl)
library(writexl)
library(haven)

library(ggplot2)
library(patchwork)
library(ggsignif)
library(gghalves)
library(ggh4x)
library(vegan)
library(compositions)

library(lme4)
library(lmerTest)
library(cowplot)
library(purrr)

setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code')
source('../R_code/myUtils.R')
set.seed(123)


load('data/beta_diversity_rawdata_forplot.RData')

bcdist_ITS = read.csv('data/Weighted_ALL_ITS.csv')
bcdist_MGS = read.csv('data/Weighted_ALL_Taxonomy.csv')
bcdist_ITS <- bcdist_ITS[!is.na(bcdist_ITS$wt1), ]
bcdist_MGS <- bcdist_MGS[!is.na(bcdist_MGS$wt2), ]


bcdist_ITS$period <- sub("^V", "T", as.character(bcdist_ITS$period))
bcdist_ITS <- bcdist_ITS %>% 
  # 1. 删除 HDP 为 NA 的行
  filter(!is.na(HDP)) %>% 
  # 2. 按规则生成 group 列
  mutate(
    group = case_when(
      HDP == 1 & preeclampsia == 1 ~ "PE",
      HDP == 1 & preeclampsia == 0 ~ "GH",
      HDP == 0                    ~ "NP",
      TRUE                        ~ NA_character_   # 兜底情况（以防出现其他取值）
    )
  )


bcdist_MGS$period <- sub("^V", "T", as.character(bcdist_MGS$period))
bcdist_MGS <- bcdist_MGS %>% 
  # 1. 删除 HDP 为 NA 的行
  filter(!is.na(HDP)) %>% 
  # 2. 按规则生成 group 列
  mutate(
    group = case_when(
      HDP == 1 & preeclampsia == 1 ~ "PE",
      HDP == 1 & preeclampsia == 0 ~ "GH",
      HDP == 0                    ~ "NP",
      TRUE                        ~ NA_character_   # 兜底情况（以防出现其他取值）
    )
  )


keep_cols_ITS <- c("id", "period", "group", all_of(fung_its_uni))

species_sample_ITS <- bcdist_ITS[, keep_cols_ITS, drop = FALSE]


keep_cols_MGS <- c("id", "period", "group", all_of(bact_mgs_s_uni))

species_sample_MGS <- bcdist_MGS[, keep_cols_MGS, drop = FALSE]


datac = species_sample_ITS %>% left_join(species_sample_MGS)
datac = datac[complete.cases(datac),]
# datac = datac %>% filter(group=='NP')


abund_tab_MGS = decostand(select(datac, all_of(bact_mgs_s_uni)), 'total', 1)

abund_tab_ITS = decostand(select(datac, all_of(fung_its_uni)), 'total', 1)

eucdist_MGS <- vegdist(abund_tab_MGS, method = "bray")

eucdist_ITS <- vegdist(abund_tab_ITS, method = "bray")


pcoa_w_unifrac_MGS = cmdscale(eucdist_MGS, k = 2, eig = TRUE)
pcoa_w_unifrac_ITS = cmdscale(eucdist_ITS, k = 2, eig = TRUE)

vpcoa = pcoa_w_unifrac_MGS

bpcoa = pcoa_w_unifrac_ITS


set.seed(9999)
# vdist = vegdist(datac %>% select(all_of(vlist)), 'bray')
# vpcoa = cmdscale(vdist, k = 2, eig = TRUE) # 只比较前2维是否相似
# bdist = vegdist(datac %>% select(all_of(blist)), 'bray')
# bpcoa = cmdscale(bdist, k = 2, eig = TRUE) # 只比较前2维是否相似


# mantelr = mantel(eucdist_MGS, eucdist_ITS, method = 'pearson', permutations = 999, parallel = 8)
# mantelr # 0.0388 0.001 
proc = protest(X = vpcoa, Y = bpcoa, permutations = 999, parallel = 8)
# summary(proc); proc



df = rbind(cbind(datac %>% select(id, period, group), type='Bacteriome', vpcoa$points),
           cbind(datac %>% select(id, period, group), type='Mycobiome', bpcoa$points))
names(df)[5:6] = c('PCoA1', 'PCoA2')
ggplot(df, aes(PCoA1, PCoA2, color = type))+
  geom_point(shape = 16, alpha = 0.8, size = 1)+
  scale_color_manual(values = c('#283618','#dda15e'))+
  scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2))+
  labs(title=expression(paste("Procrustes ",italic("R"), "=0.10, ", italic("P"), "< 0.001 (ALL)")))+
  theme_classic(base_line_size = 0.3)+
  theme(axis.text = element_text(color='black', size=12),
        plot.title.position = "plot",
        plot.title = element_text(size = 11, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=11),
        legend.key.spacing = unit(1.2,'mm'),
        legend.key.size = unit(3,'mm'),
        legend.position = c(1,0), legend.justification = c(1, 0),
        legend.background = element_rect(fill='transparent'),
        plot.margin = ggplot2::margin(1,3,1,1,'mm'))
ggsave('Anti_graph/fig2e_pcoa_intera_ALL.pdf', width = 3, height = 2.8)


