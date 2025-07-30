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
library(purrr)

setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code')
source('../R_code/myUtils.R')
set.seed(123)



load('../data/beta_diversity_rawdata_forplot.RData')

## beta diversity plots----

ptColors2 = c('#7a94ba','#efb0b3')
ptColors3 = c('#7a94ba','#f8c672','#B5665D')

mytheme <- theme_bw()+
  theme(plot.background = element_blank(),  # 移除图表背景
        panel.background = element_blank(),  # 移除绘图区域背景
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 18, color = "black", vjust = 0.5), 
        axis.text.y = element_text(size = 18, color = "black", hjust = 0.5),
        strip.text = element_text(size = 14, face = "bold"),
        axis.ticks.length = unit(1.5, 'mm'),
        title = element_text(size = 8.7),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = c(0.95, 0.11),
        legend.background = element_rect(fill = NA),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),  # 调整整体边距
        panel.spacing = unit(0.1, "cm"),  # 减小面板之间的间距
        axis.line = element_line(colour = "black", size = 0.5),  # 添加坐标轴线
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5) 
  )


# 时期 term
# 输入残差矩阵 in_bcdist
# 菌群名称 in_trans
# 保存路径files_path

beta_plot <- function(term, bcdist, trans, files_path){
  # term = 'V1'
  # bcdist = bact_mgs_s_dist
  # trans = bact_mgs_s_uni
  
  bcdist <- bcdist %>% filter(period == term)
  eucdist <- bcdist %>% select(all_of(trans)) %>% vegdist('euclidean')   
  pcoa <- cmdscale(eucdist, k = 3, eig = TRUE) #只存前3维
  contr <- pcoa$eig/sum(pcoa$eig) * 100
  colnames(pcoa$points) = c('PCOA1', 'PCOA2', 'PCOA3')
  pbeta <- cbind(bcdist, pcoa$points[,1:3])
  
  test <- adonis2(eucdist ~ HDP, data = bcdist, permutations = 999, parallel = 8)
  adonis <- paste(str_replace(term, "V", "T"),"\nR2 = ", round(test$R2[1],3), "\nP = ", round(test$`Pr(>F)`[1],3))
  
  pbeta$HDP_label <- factor(pbeta$HDP, levels = c(0, 1), labels = c("Non-HDP", "HDP"))
  
  # main plot
  main <- 
    ggplot(pbeta, aes(x=PCOA1, y=PCOA2))+ 
    geom_point(aes(color = HDP_label), size=1.5, alpha=0.8, shape = 16)+
    stat_ellipse(aes(color = HDP_label, fill=HDP_label), geom = "polygon",level = 0.95,
                 linetype = 1, linewidth=0.8, alpha = 0.1)+
    scale_colour_manual(values=ptColors2)+
    scale_x_continuous(expand = c(0.02, 0.02)) +  # 添加这行
    scale_y_continuous(expand = c(0.02, 0.02)) +  # 添加这行
    scale_fill_manual(values=ptColors2)+
    labs(x=sprintf('PCoA1 (%.1f%%)', contr[1]), 
         y=sprintf('PCoA2 (%.1f%%)', contr[2]),
         color = 'Group', fill = 'Group')+
    mytheme+
    # coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10))+     theme(legend.position = c(0.8, 0.2))
    theme(legend.position = "none")

  box1 <- ggplot(pbeta, 
                 aes(x = HDP_label, y = PCOA1, fill = HDP_label)) +
    geom_boxplot(show.legend = FALSE) + 
    stat_boxplot(geom = "errorbar", width = 0.1, size = 0.5) +
    #geom_jitter(show.legend = FALSE) + 
    scale_fill_manual(values = ptColors2) +
    coord_flip() + 
    theme_bw() +
    theme(panel.grid=element_blank(), 
          axis.title = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(color='black'),
          axis.text.x = element_blank(),
          axis.text.y = element_text(colour = 'black', size = 14))
  
  box2 <- ggplot(pbeta, 
                 aes(x = HDP_label, y = PCOA2, fill = HDP_label)) +
    geom_boxplot(show.legend = FALSE) + 
    stat_boxplot(geom = "errorbar", width = 0.1, size = 0.5) +
    #geom_jitter(show.legend = FALSE) + 
    scale_fill_manual(values = ptColors2) +
    theme_bw() +
    theme(panel.grid=element_blank(), 
          axis.title = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(color='black'),
          axis.text.y = element_blank(),
          axis.text.x =element_text(colour ='black', size = 14,  angle = 45, vjust = 1, hjust = 1))
  
  box3 <- ggplot(pbeta, aes(PCOA1, PCOA2)) +
    annotate("text", x = -0.5, y = 0.2, 
             label = adonis,  # 直接使用原始文本
             size = 5, 
             fontface = "bold") +  # 使用 fontface 设置粗体
    theme_bw() +
    xlab("") + ylab("") +
    theme(panel.grid=element_blank(), 
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())
  
  title <- ggplot() +
    annotate("text", x = 0.2, y = 0.8, label = "Metagenomic species", size = 8, fontface = "bold") +
    theme_void()
  
  p <- title + (box1 + box3 + main + box2 +
                  plot_layout(heights = c(1,4), widths = c(4,1), ncol = 2, nrow = 2)) +
    plot_layout(heights = c(0.8, 10), ncol = 1)  # 调整标题和主图的比例
  
  ggsave(plot = p, width = 8, height = 8, device = "pdf", file = files_path)
  
  
}

beta_plot('V3', bact_mgs_s_dist, bact_mgs_s_uni, 'Anti_graph/Beta/Beta_MGS_S_V3.pdf') #switch to bact_mgs_g_dist bact_16s_dist fung_its_dist





