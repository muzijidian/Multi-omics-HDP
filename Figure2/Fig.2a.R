rm(list = ls())
library(pheatmap)
library(openxlsx)
library(tidyr)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(stringr)
library(readxl)
library(stringi)
library(xml2)

setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code')
source('../R_code/myUtils.R')
set.seed(123)


library(haven)




#########################
HDP_type = "health"
type_list = c('Metabolome', 'ITS', "16S", 'Metagenomic genus', 'Metagenomic species')
if(HDP_type == "diease"){
  df = read.xlsx('../data/omics_data_timevar_icc_HDP.xlsx')
}else{
  df = read.xlsx('../data/omics_data_timevar_icc_NonHDP.xlsx')
}


df = df %>%
  mutate(data_type = case_when(
    data_type == "Taxonomy_s" ~ "Metagenomic species",
    data_type == "Taxonomy_g" ~ "Metagenomic genus",
    TRUE ~ data_type
  )) %>%
  mutate(class = case_match(
    data_type,
    'Metabolome' ~ 'Metabolome',
    c('16S', 'ITS') ~ '16S and ITS',
    c('Metagenomic genus', 'Metagenomic species') ~ 'Metagenomic taxonomy'
  ))
df$data_type = factor(df$data_type, levels = type_list)
df$class = factor(df$class, levels = c('Metabolome','16S and ITS','Metagenomic taxonomy'))

max_row <- df %>%
  filter(data_type == "Metagenomic genus") %>%
  slice_max(order_by = ICC, n = 1)

# 查看结果
print(max_row)

# type_col = c("#839B97", "#FF8882", "#6E7582", "#B6A39E", "#4A7B9D")
type_col = c("#472D7B3a", "#ebdcb2", "#95C6C6", "#FFB2B9", "#7EABCA")

p <- ggplot(df, aes(x=time_var, y=ICC, color=data_type, fill=data_type))+ 
  facet_wrap(~class, scales = 'free')+
  geom_point(size=3, alpha=0.7, shape = 16)+
  scale_color_manual(values=type_col)+
  labs(x = 'Variance explained by gestational weeks')+
  labs(y = 'Variance attributed to participants')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(fill = "#d3dce3", color = 'black'),## #d3dce3 f4e1da
        strip.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, color = "black", face = "bold"),
        axis.text = element_text(size = 20, color = 'black'),
        legend.text = element_text(size = 18, face = "bold"), 
        legend.key.size = unit(1.2, "cm"),
        legend.title = element_blank())












