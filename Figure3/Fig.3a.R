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
# library(ComplexHeatmap)
# library(circlize)

setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code')
source('../R_code/myUtils.R')
set.seed(123)




########################Metagenome HDP######################

data_type = "Taxonomy"
level = 's'
transformer = "c"
wt = "wt2"

state = "all"



file_name = paste0('../data/HDP_ALL_', data_type, '_', level, '_', transformer, '_', wt, '_AbundInfo0.25_Info_no_antibiotic_as_cov.xlsx')

HDP_T1 = read.xlsx(file_name, sheet = 1)
HDP_T2 = read.xlsx(file_name, sheet = 2)
HDP_T3 = read.xlsx(file_name, sheet = 3)
PE_T1 = read.xlsx(file_name, sheet = 4)
PE_T2 = read.xlsx(file_name, sheet = 5)
PE_T3 = read.xlsx(file_name, sheet = 6)
GH_T1 = read.xlsx(file_name, sheet = 7)
GH_T2 = read.xlsx(file_name, sheet = 8)
GH_T3 = read.xlsx(file_name, sheet = 9)

HDP_T1$source <- "HDP_T1"
HDP_T2$source <- "HDP_T2"
HDP_T3$source <- "HDP_T3"
PE_T1$source <- "PE_T1"
PE_T2$source <- "PE_T2"
PE_T3$source <- "PE_T3"
GH_T1$source <- "GH_T1"
GH_T2$source <- "GH_T2"
GH_T3$source <- "GH_T3"



if(state=="T1"){
  combined_data <- bind_rows(HDP_T1, PE_T1, GH_T1)
}
if(state=="T23"){
  combined_data <- bind_rows(HDP_T2, PE_T2, GH_T2, HDP_T3,  PE_T3, GH_T3)
}
if(state=="all"){
  combined_data <- bind_rows(HDP_T1, PE_T1, GH_T1, HDP_T2, PE_T2, GH_T2, HDP_T3, PE_T3, GH_T3)
}




file_name2 = paste0('../data/HDP_ALL_', data_type, '_g_', transformer, '_', wt, '_AbundInfo0.25_Info_no_antibiotic_as_cov.xlsx')
print(file_name2)


HDP_T1_2 = read.xlsx(file_name2, sheet = 1)
HDP_T1_2$source <- "HDP_T1"
HDP_T2_2 = read.xlsx(file_name2, sheet = 2)
HDP_T2_2$source <- "HDP_T2"
HDP_T3_2 = read.xlsx(file_name2, sheet = 3)
HDP_T3_2$source <- "HDP_T3"
PE_T1_2 = read.xlsx(file_name2, sheet = 4)
PE_T1_2$source <- "PE_T1"
PE_T2_2 = read.xlsx(file_name2, sheet = 5)
PE_T2_2$source <- "PE_T2"
PE_T3_2 = read.xlsx(file_name2, sheet = 6)
PE_T3_2$source <- "PE_T3"
GH_T1_2 = read.xlsx(file_name2, sheet = 7)
GH_T1_2$source <- "GH_T1"
GH_T2_2 = read.xlsx(file_name2, sheet = 8)
GH_T2_2$source <- "GH_T2"
GH_T3_2 = read.xlsx(file_name2, sheet = 9)
GH_T3_2$source <- "GH_T3"


combined_data2 = bind_rows(HDP_T1_2, PE_T1_2, GH_T1_2, HDP_T2_2, PE_T2_2, GH_T2_2, HDP_T3_2, PE_T3_2, GH_T3_2)
combined_data2 <- combined_data2[combined_data2$meta_name != "g93", ] # since it is unclassified
rownames(combined_data2) <- NULL


combined_data <- bind_rows(combined_data, combined_data2)


g_data <- combined_data %>% filter(str_starts(name, "g__"))
s_data <- combined_data %>% filter(str_starts(name, "s__"))




g_priority_names <- g_data %>%
  group_by(name) %>%
  filter(any(source %in% c("HDP_T1", "PE_T1", "GH_T1") & pv_adj_FDR < 0.1)) %>%
  ungroup() %>%
  distinct(name) %>%
  arrange(name) %>%
  pull(name)


g_priority_df <- g_data %>%
  filter(name %in% g_priority_names) %>%
  distinct(name, .keep_all = TRUE) %>%  
  select(name, meta_name)


g_remaining_names <- g_data %>%
  filter(!(name %in% g_priority_names)) %>%
  distinct(name) %>%
  arrange(name) %>%
  pull(name)


s_priority_names <- s_data %>%
  group_by(name) %>%
  filter(any(source %in% c("HDP_T1", "PE_T1", "GH_T1") & pv_adj_FDR < 0.1)) %>%
  ungroup() %>%
  distinct(name) %>%
  arrange(name) %>%
  pull(name)

s_priority_df <- s_data %>%
  filter(name %in% s_priority_names) %>%
  distinct(name, .keep_all = TRUE) %>%  
  select(name, meta_name)


s_remaining_names <- s_data %>%
  filter(!(name %in% s_priority_names)) %>%
  distinct(name) %>%
  arrange(name) %>%
  pull(name)


final_row_names <- c(g_priority_names, g_remaining_names, s_priority_names, s_remaining_names)


unique_meta_names <- unique(combined_data$meta_name)
print(unique_meta_names)


pv_data <- combined_data %>%
  select(name, source, pv_adj_FDR) %>%
  spread(key = source, value = pv_adj_FDR)
row_names <- pv_data$name 
pv_data$name <- NULL
rownames(pv_data) <- row_names

orig_pv_data <- combined_data %>%
  select(name, source, pv) %>%
  spread(key = source, value = pv)
row_names <- orig_pv_data$name 
orig_pv_data$name <- NULL
rownames(orig_pv_data) <- row_names



beta_data <- combined_data %>%
  select(name, source, beta, phylum) %>%
  pivot_wider(names_from = source, values_from = beta)



beta_data <- beta_data %>%
  mutate(name = factor(name, levels = final_row_names)) %>%
  arrange(name) 
beta_data$name <- NULL
rownames(beta_data) <- final_row_names


beta_data_orig <- beta_data%>%select(-phylum)
rownames(beta_data_orig) <- final_row_names

row_annotation <- data.frame(phylum = beta_data$phylum)
rownames(row_annotation) <- rownames(beta_data_orig)


unique_phylum <- unique(beta_data$phylum)
print(unique_phylum)


if (level == 'g'){
  
  if(state=="T1"){
    unique_phylum <- unique(row_annotation$phylum)
    phylum_colors <- setNames(colorRampPalette(brewer.pal(3, "Set3"))(length(unique_phylum)), unique_phylum)
  } else if(state=="T23"){
    unique_phylum <- c("Bacteroidetes", "Firmicutes", "Proteobacteria", 
                       "Actinobacteria") 
    phylum_colors <- setNames(colorRampPalette(brewer.pal(4, "Set3"))(length(unique_phylum)), unique_phylum)
  } else if(state=="all"){
    unique_phylum <- c("Bacteroidetes", "Firmicutes", "Proteobacteria", 
                       "Actinobacteria") 
    phylum_colors <- setNames(colorRampPalette(brewer.pal(4, "Set3"))(length(unique_phylum)), unique_phylum)
  }
}



if (level == 's'){
  
  if(state=="T1"){
    unique_phylum <- unique(row_annotation$phylum)
    phylum_colors <- setNames(colorRampPalette(brewer.pal(3, "Set3"))(length(unique_phylum)), unique_phylum)
  } else if(state=="T23"){
    unique_phylum <- c("Bacteroidetes", "Firmicutes", "Proteobacteria", 
                       "Actinobacteria") 
    phylum_colors <- setNames(colorRampPalette(brewer.pal(4, "Set3"))(length(unique_phylum)), unique_phylum)
  } else if(state=="all"){
    unique_phylum <- c("Bacteroidetes", "Firmicutes", "Proteobacteria", 
                       "Actinobacteria") 
    phylum_colors <- setNames(colorRampPalette(brewer.pal(4, "Set3"))(length(unique_phylum)), unique_phylum)
  }
}


annotation_colors <- list(phylum = phylum_colors)



heatmap_result <- pheatmap(beta_data_orig, 
                           cluster_cols = FALSE, 
                           annotation_row = row_annotation, 
                           annotation_colors = annotation_colors)

column_order <- colnames(beta_data_orig)

pv_data <- pv_data[final_row_names, column_order, drop = FALSE]
rownames(pv_data) <- final_row_names


orig_pv_data <- orig_pv_data[final_row_names, column_order, drop = FALSE]

rownames(orig_pv_data) <- final_row_names


beta_data_orig[orig_pv_data >= 0.05] <- NA
rownames(beta_data_orig) <- rownames(orig_pv_data)


significance_matrix <- ifelse(pv_data < 0.01, "**",
                              ifelse(pv_data < 0.05, "*",
                                     ifelse(pv_data < 0.1, "+", "")))

n_cols <- ncol(beta_data_orig)
gaps <- seq(3, n_cols - 1, by = 3)

row_labels_colors <- ifelse(startsWith(rownames(beta_data_orig), "g__"), "#edf5f5", "#f4f2ec")

library(grid)

if(level=='s'){
  beta_data_orig <- apply(beta_data_orig, c(1, 2), function(x) pmax(pmin(x, 1.), -1.))

  breaks <- seq(-1., 1., length.out = 101)

  heatmap <- pheatmap(
    beta_data_orig,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = row_annotation,
    annotation_colors = annotation_colors,
    display_numbers = significance_matrix,
    number_color = "blue",
    gaps_col = gaps,
    angle_col = "315",
    na_col = "white",
    background = "transparent",
    breaks = breaks,
    width = 7.1, height = 11.7,
    legend_position = "bottomright",
    legend = TRUE,
    fontsize_row = 11, 
    fontsize_col = 11
  )
}










############HDP 16S ITS##################################



data_type = "16S"
transformer = "c"
wt = "wt1"
state = "all"




HDP_T1 = read.xlsx(paste0('../data/HDP_ALL_', data_type, '_', transformer, '_', wt, '_AbundInfo0.25_Info_no_antibiotic_as_cov.xlsx'), sheet = 1)
HDP_T2 = read.xlsx(paste0('../data/HDP_ALL_', data_type, '_', transformer, '_', wt, '_AbundInfo0.25_Info_no_antibiotic_as_cov.xlsx'), sheet = 2)
HDP_T3 = read.xlsx(paste0('../data/HDP_ALL_', data_type, '_', transformer, '_', wt, '_AbundInfo0.25_Info_no_antibiotic_as_cov.xlsx'), sheet = 3)
PE_T1 = read.xlsx(paste0('../data/HDP_ALL_', data_type, '_', transformer, '_', wt, '_AbundInfo0.25_Info_no_antibiotic_as_cov.xlsx'), sheet = 4)
PE_T2 = read.xlsx(paste0('../data/HDP_ALL_', data_type, '_', transformer, '_', wt, '_AbundInfo0.25_Info_no_antibiotic_as_cov.xlsx'),, sheet = 5)
PE_T3 = read.xlsx(paste0('../data/HDP_ALL_', data_type, '_', transformer, '_', wt, '_AbundInfo0.25_Info_no_antibiotic_as_cov.xlsx'), sheet = 6)
GH_T1 = read.xlsx(paste0('../data/HDP_ALL_', data_type, '_', transformer, '_', wt, '_AbundInfo0.25_Info_no_antibiotic_as_cov.xlsx'), sheet = 7)
GH_T2 = read.xlsx(paste0('../data/HDP_ALL_', data_type, '_', transformer, '_', wt, '_AbundInfo0.25_Info_no_antibiotic_as_cov.xlsx'), sheet = 8)
GH_T3 = read.xlsx(paste0('../data/HDP_ALL_', data_type, '_', transformer, '_', wt, '_AbundInfo0.25_Info_no_antibiotic_as_cov.xlsx'), sheet = 9)

HDP_T1$source <- "HDP_T1"
HDP_T2$source <- "HDP_T2"
HDP_T3$source <- "HDP_T3"
PE_T1$source <- "PE_T1"
PE_T2$source <- "PE_T2"
PE_T3$source <- "PE_T3"
GH_T1$source <- "GH_T1"
GH_T2$source <- "GH_T2"
GH_T3$source <- "GH_T3"



if(state=="T1"){
  combined_data <- bind_rows(HDP_T1, PE_T1, GH_T1)
}
if(state=="T23"){
  combined_data <- bind_rows(HDP_T2, PE_T2, GH_T2, HDP_T3,  PE_T3, GH_T3)
}
if(state=="all"){
  combined_data <- bind_rows(HDP_T1, PE_T1, GH_T1, HDP_T2, PE_T2, GH_T2, HDP_T3, PE_T3, GH_T3)
}

combined_data <- combined_data %>%
  mutate(data_source = '16S')




file_name2 = paste0('../data/HDP_ALL_ITS_c_wt1_AbundInfo0.25_Info_no_antibiotic_as_cov.xlsx')

HDP_T1_2 = read.xlsx(file_name2, sheet = 1)
HDP_T1_2$source <- "HDP_T1"
HDP_T2_2 = read.xlsx(file_name2, sheet = 2)
HDP_T2_2$source <- "HDP_T2"
HDP_T3_2 = read.xlsx(file_name2, sheet = 3)
HDP_T3_2$source <- "HDP_T3"
PE_T1_2 = read.xlsx(file_name2, sheet = 4)
PE_T1_2$source <- "PE_T1"
PE_T2_2 = read.xlsx(file_name2, sheet = 5)
PE_T2_2$source <- "PE_T2"
PE_T3_2 = read.xlsx(file_name2, sheet = 6)
PE_T3_2$source <- "PE_T3"
GH_T1_2 = read.xlsx(file_name2, sheet = 7)
GH_T1_2$source <- "GH_T1"
GH_T2_2 = read.xlsx(file_name2, sheet = 8)
GH_T2_2$source <- "GH_T2"
GH_T3_2 = read.xlsx(file_name2, sheet = 9)
GH_T3_2$source <- "GH_T3"


combined_data2 = bind_rows(HDP_T1_2, PE_T1_2, GH_T1_2, HDP_T2_2, PE_T2_2, GH_T2_2, HDP_T3_2, PE_T3_2, GH_T3_2)
combined_data2 <- combined_data2 %>%
  mutate(data_source = 'ITS')



combined_data <- bind_rows(combined_data2, combined_data)


combined_data <- combined_data %>%
  mutate(name = ifelse(name == "g__uncultured", meta_name, name))
combined_data <- combined_data %>%
  mutate(beta = pmin(pmax(beta, -3), 3))

combined_data <- combined_data[!grepl('g[0-9]', combined_data$name), ]

g_data <- combined_data %>% filter(data_source == "ITS")
s_data <- combined_data %>% filter(data_source == "16S")


g_priority_names <- g_data %>%
  group_by(name) %>%
  filter(any(source %in% c("HDP_T2", "PE_T2", "GH_T2") & pv_adj_FDR < 0.1)) %>%
  ungroup() %>%
  distinct(name) %>%
  arrange(name) %>%  
  pull(name)

g_priority_df <- g_data %>%
  filter(name %in% g_priority_names) %>%
  distinct(name, .keep_all = TRUE) %>%  
  select(name, meta_name)


g_remaining_names <- g_data %>%
  filter(!(name %in% g_priority_names)) %>%
  distinct(name) %>%
  arrange(name) %>%  
  pull(name)


s_priority_names <- s_data %>%
  group_by(name) %>%
  filter(any(source %in% c("HDP_T1", "PE_T1", "GH_T1") & pv_adj_FDR < 0.1)) %>%
  ungroup() %>%
  distinct(name) %>%
  arrange(name) %>% 
  pull(name)

s_priority_df <- s_data %>%
  filter(name %in% s_priority_names) %>%
  distinct(name, .keep_all = TRUE) %>%  
  select(name, meta_name)



s_remaining_names <- s_data %>%
  filter(!(name %in% s_priority_names)) %>%
  distinct(name) %>%
  arrange(name) %>%  
  pull(name)


final_row_names <- c(g_priority_names, g_remaining_names, s_priority_names, s_remaining_names)



unique_meta_names <- unique(combined_data$meta_name)
print(unique_meta_names)




pv_data <- combined_data %>%
  select(name, source, pv_adj_FDR) %>%
  spread(key = source, value = pv_adj_FDR)
row_names <- pv_data$name 
pv_data$name <- NULL
rownames(pv_data) <- row_names

orig_pv_data <- combined_data %>%
  select(name, source, pv) %>%
  spread(key = source, value = pv)
row_names <- orig_pv_data$name 
orig_pv_data$name <- NULL
rownames(orig_pv_data) <- row_names



beta_data <- combined_data %>%
  select(name, source, beta, phylum) %>%
  pivot_wider(names_from = source, values_from = beta)



beta_data <- beta_data %>%
  mutate(name = factor(name, levels = final_row_names)) %>%
  arrange(name) 
beta_data$name <- NULL
rownames(beta_data) <- final_row_names


beta_data_orig <- beta_data%>%select(-phylum)
rownames(beta_data_orig) <- final_row_names

row_annotation <- data.frame(phylum = beta_data$phylum)
rownames(row_annotation) <- rownames(beta_data_orig)

unique_phylum <- unique(beta_data$phylum)
print(unique_phylum)

if(data_type=="16S"){

if(state=="T1"){
unique_phylum <- unique(row_annotation$phylum)
phylum_colors <- setNames(colorRampPalette(brewer.pal(3, "Set3"))(length(unique_phylum)), unique_phylum)
}

if(state=="T23"){
unique_phylum <- c("Bacteroidota", "Firmicutes", "Proteobacteria", 
                   "Actinobacteriota", "Cyanobacteria") 
phylum_colors <- setNames(colorRampPalette(brewer.pal(5, "Set3"))(length(unique_phylum)), unique_phylum)
}

if(state=="all"){
  unique_phylum <- c("Bacteroidota", "Firmicutes", "Proteobacteria", "Ascomycota", "Basidiomycota",
                     "Actinobacteriota", "Cyanobacteria", "Verrucomicrobiota", "Desulfobacterota") ##"Ascomycota", "Basidiomycota", "Mucoromycota"
  phylum_colors <- setNames(colorRampPalette(brewer.pal(9, "Set3"))(length(unique_phylum)), unique_phylum)
}
}

if(data_type=="ITS"){
  
  if(state=="T1"){
    unique_phylum <- unique(row_annotation$phylum)
    phylum_colors <- setNames(colorRampPalette(brewer.pal(3, "Set3"))(length(unique_phylum)), unique_phylum)
  }
  
  if(state=="T23"){
    unique_phylum <- unique(row_annotation$phylum)
    phylum_colors <- setNames(colorRampPalette(brewer.pal(3, "Set3"))(length(unique_phylum)), unique_phylum)
  }
  
  if(state=="all"){
    unique_phylum <- unique(row_annotation$phylum)
    phylum_colors <- setNames(colorRampPalette(brewer.pal(3, "Set3"))(length(unique_phylum)), unique_phylum)
  }
}

annotation_colors <- list(phylum = phylum_colors)


heatmap_result <- pheatmap(beta_data_orig, 
                           cluster_cols = FALSE, 
                           annotation_row = row_annotation, 
                           annotation_colors = annotation_colors)


print(final_row_names)
column_order <- colnames(beta_data_orig)

pv_data <- pv_data[final_row_names, column_order, drop = FALSE]
rownames(pv_data) <- final_row_names


orig_pv_data <- orig_pv_data[final_row_names, column_order, drop = FALSE]


rownames(orig_pv_data) <- final_row_names


beta_data_orig[orig_pv_data >= 0.05] <- NA
rownames(beta_data_orig) <- rownames(orig_pv_data)


significance_matrix <- ifelse(pv_data < 0.01, "**",
                              ifelse(pv_data < 0.05, "*",
                                     ifelse(pv_data < 0.1, "+", "")))

n_cols <- ncol(beta_data_orig)

gaps <- seq(3, n_cols - 1, by = 3)


if(data_type=='16S'){
  
  beta_data_orig <- apply(beta_data_orig, c(1, 2), function(x) pmax(pmin(x, 1), -1))

  breaks <- seq(-1, 1, length.out = 101)

  pheatmap(beta_data_orig, 
           cluster_rows = FALSE, 
           cluster_cols = FALSE, 
           annotation_row = row_annotation, 
           annotation_colors = annotation_colors, 
           display_numbers = significance_matrix,  # 添加显著性标注
           number_color = "blue",
           gaps_col = gaps,
           angle_col = "315",
           na_col = "white",
           breaks = breaks,
           # width = 6.3, height =3.5, # for T1
           width = 7.83, height =9.754, # for T23
           legend_position = "bottomright",
           legend = TRUE)
}


