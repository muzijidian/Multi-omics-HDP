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




#############################Shannon################################

load('../data/alpha_diversity_raw.RData')

data_16S <-  alpha_fung_its  # switch to alpha_bact_16s  alpha_bact_mgs_g alpha_bact_mgs_s
data_16S$source <- "ITS" 

data_16S$group <- ifelse(data_16S$HDP == 1, "HDP", "Non-HDP")

filtered_data <- data_16S[data_16S$index == "Shannon", ] # switch to Simpson/Chao1
filtered_data <- filtered_data %>%
  mutate(period = recode(period, "V1" = "T1", "V2" = "T2", "V3" = "T3"))

names(filtered_data)[names(filtered_data) == "alpha"] <- "Shannon"


calculate_trend_pv <- function(alpha_div, group, metric){
  df_trend = alpha_div %>% mutate(period = as.numeric(str_remove(period, 'T')))
  rule = as.formula(paste0(metric, '~(1|id)+period'))
  fit = lmer(rule, data = df_trend[df_trend[,group]==1,], REML = F)
  round(summary(fit)$coefficients['period', ], 4) %>% print()
  fit = lmer(rule, data = df_trend[df_trend[,group]==0,], REML = F)
  round(summary(fit)$coefficients['period', ], 4) %>% print()
  fit1 = lmer(Shannon~(1|id)+HDP*period, data=df_trend, REML=F)
  round(summary(fit1)$coefficients[paste0(group,':period'), ], 4) %>% print()
}



calculate_trend_pv(filtered_data, 'HDP', 'Shannon')

col_group = c('HDP'='#cc8082','Non-HDP'='#5c6d8e')

p_shannon = ggplot(filtered_data, aes(x = period, y = Shannon, fill = group, color=group))+
  facet_wrap(~"Shannon")+
  geom_violin(width = 0.3, alpha = 0.2, position = position_dodge(0.4))+
  geom_boxplot(size= 0.8, width = 0.3, alpha = 0.1, outlier.shape = NA, 
               position = position_dodge(0.4))+
  # geom_point(alpha=0.5, size=1,
  #            position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0,
  #                                            dodge.width = 0.4))+
  stat_summary(fun = median, geom = "line", aes(group = group), 
               position = position_dodge(0.4))+
  scale_fill_manual(values = col_group)+
  scale_color_manual(values = col_group)+
  # scale_y_continuous(expand = c(0,0.1))+
  scale_y_continuous(expand = c(0,0.1), limits = c(NA, 4.3))+
  annotate('text', x = 0.55, y = 4.16, size = 3.35, hjust = 0,
           label=expression(paste(italic('P'),'-interaction = 0.015')))+
  annotate('text', x = 0.55, y = 3.875, size = 3.35, hjust = 0,
           label=expression(paste(italic('P'),'-trend = 0.006 for HDP')))+
  annotate('text', x = 0.55, y = 3.60, size = 3.35, hjust = 0,
           label=expression(paste(italic('P'),'-trend = 0.446 for Non-HDP')))+
  stat_compare_means(aes(group = group), method = "wilcox.test",
                     label = "p.signif", show.legend = F, label.y = 3.2)+
  theme_bw()+
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 10.2, color = "black"), 
    axis.text.y = element_text(size = 9.2, color = "grey10"),
    strip.text = element_text(size = 11, color = 'black'),
    strip.background = element_rect(fill = '#c0a8b6'),
    legend.title = element_blank(),
    panel.grid= element_line(linetype = 2),
    legend.position = "none" 
  )




