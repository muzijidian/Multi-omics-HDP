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
level = 'g'
transformer = "c"
wt = "wt2"

state = "all"



file_name = paste0('../Anti_data/preprocessed/HDP_ALL_', data_type, '_', level, '_', transformer, '_', wt, '_AbundInfo0.25_Info_no_antibiotic_as_cov_class3.xlsx')

# HDP_T1 = read.xlsx(file_name, sheet = 1)
# HDP_T2 = read.xlsx(file_name, sheet = 2)
# HDP_T3 = read.xlsx(file_name, sheet = 3)
# HDP_Lmer = read.xlsx(file_name, sheet = 10)
PE_T1 = read.xlsx(file_name, sheet = 4)
PE_T2 = read.xlsx(file_name, sheet = 5)
PE_T3 = read.xlsx(file_name, sheet = 6)
# PE_Lmer = read.xlsx(file_name, sheet = 11)
GH_T1 = read.xlsx(file_name, sheet = 7)
GH_T2 = read.xlsx(file_name, sheet = 8)
GH_T3 = read.xlsx(file_name, sheet = 9)
# GH_Lmer = read.xlsx(file_name, sheet = 12)
PEvsGH_T1 = read.xlsx(file_name, sheet = 10)
PEvsGH_T2 = read.xlsx(file_name, sheet = 11)
PEvsGH_T3 = read.xlsx(file_name, sheet = 12)

PE_Pooled = read.xlsx(file_name, sheet = 13)
GH_Pooled = read.xlsx(file_name, sheet = 14)
PEvsGH_Pooled = read.xlsx(file_name, sheet = 15)

# HDP_T1$source <- "HDP_T1"
# HDP_T2$source <- "HDP_T2"
# HDP_T3$source <- "HDP_T3"
PE_T1$source <- "PE_T1"
PE_T2$source <- "PE_T2"
PE_T3$source <- "PE_T3"
GH_T1$source <- "GH_T1"
GH_T2$source <- "GH_T2"
GH_T3$source <- "GH_T3"
PEvsGH_T1$source <- "PE vs. GH (T1)"
PEvsGH_T2$source <- "PE vs. GH (T2)"
PEvsGH_T3$source <- "PE vs. GH (T3)"

PE_Pooled$source <- "Pooled (PE)"
GH_Pooled$source <- "Pooled (GH)"
PEvsGH_Pooled$source <- "Pooled (PE vs. GH)"


if(state=="T1"){
  combined_data <- bind_rows(HDP_T1, PE_T1, GH_T1)
}
if(state=="T23"){
  combined_data <- bind_rows(HDP_T2, PE_T2, GH_T2, HDP_T3,  PE_T3, GH_T3)
}
if(state=="all"){
  combined_data <- bind_rows(PE_T1, GH_T1, PEvsGH_T1, PE_T2, GH_T2, PEvsGH_T2, PE_T3, GH_T3, PEvsGH_T3, PE_Pooled, GH_Pooled, PEvsGH_Pooled)
}

# 将 combined_data 分成以 g__ 和 s__ 开头的数据
g_data <- combined_data %>% filter(str_starts(name, "g__"))
s_data <- combined_data %>% filter(str_starts(name, "s__"))
# library(writexl)
# write_xlsx(g_data, "Anti_results/Taxonomy_g_sig_info.xlsx")
# write_xlsx(s_data, "Anti_results/Taxonomy_s_sig_info.xlsx")


# 找到以 g__ 开头并且在 HDP_T1, PE_T1, GH_T1 中任意一列 p 值小于 0.05 的 name 顺序
g_priority_names <- g_data %>%
  group_by(name) %>%
  filter(any(source %in% c("PE_T1", "GH_T1", "PEvsGH_T1") & pv_adj_FDR < 0.1)) %>%
  ungroup() %>%
  distinct(name) %>%
  arrange(name) %>%
  pull(name)


g_priority_df <- g_data %>%
  filter(name %in% g_priority_names) %>%
  distinct(name, .keep_all = TRUE) %>%  
  select(name, meta_name)


# 找到剩余的 g__ 开头的 name 顺序
g_remaining_names <- g_data %>%
  filter(!(name %in% g_priority_names)) %>%
  distinct(name) %>%
  arrange(name) %>%
  pull(name)

# 找到以 s__ 开头并且在 HDP_T1, PE_T1, GH_T1 中任意一列 p 值小于 0.05 的 name 顺序
s_priority_names <- s_data %>%
  group_by(name) %>%
  filter(any(source %in% c("PE_T1", "GH_T1", "PEvsGH_T1") & pv_adj_FDR < 0.1)) %>%
  ungroup() %>%
  distinct(name) %>%
  arrange(name) %>%
  pull(name)

s_priority_df <- s_data %>%
  filter(name %in% s_priority_names) %>%
  distinct(name, .keep_all = TRUE) %>%  
  select(name, meta_name)


write_xlsx(s_priority_df, "Anti_results/Taxonomy_s_sig.xlsx")

# 找到剩余的 s__ 开头的 name 顺序
s_remaining_names <- s_data %>%
  filter(!(name %in% s_priority_names)) %>%
  distinct(name) %>%
  arrange(name) %>%
  pull(name)

# 最终的 name 顺序
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
save_file = paste0('Anti_graph/HDP_', data_type, '_', transformer, '_', wt, '_heatmap0.25_combine_as_cov_class3.pdf')

row_labels_colors <- ifelse(startsWith(rownames(beta_data_orig), "g__"), "#edf5f5", "#f4f2ec")

library(grid)


if(level=='s'){
  beta_data_orig <- apply(beta_data_orig, c(1, 2), function(x) pmax(pmin(x, 1.), -1.))
  #
  # breaks <- seq(-0.4, 0.4, length.out = 101) #for ITS
  breaks <- seq(-1, 1, length.out = 101)

  # heatmap <- pheatmap(
  #   beta_data_orig,
  #   cluster_rows = FALSE,
  #   cluster_cols = FALSE,
  #   annotation_row = row_annotation,
  #   annotation_colors = annotation_colors,
  #   display_numbers = significance_matrix,
  #   fontsize_number = 13,
  #   number_color = "blue",
  #   gaps_col = gaps,
  #   angle_col = "315",
  #   na_col = "white",
  #   background = "transparent",
  #   breaks = breaks,
  #   filename = save_file,
  #   width = 8.2, height = 7.6,
  #   legend_position = "bottomright",
  #   legend = TRUE,
  #   fontsize_row = 11,  # 设置行标签(纵轴)字体大小
  #   fontsize_col = 11
  # )
  
  clean_row_labels <- gsub("^s__", "", rownames(beta_data_orig))
  
  # --- 改动 2: 生成热图对象，但先不保存 (filename = NA) ---
  heatmap <- pheatmap(
    beta_data_orig,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = row_annotation,
    annotation_colors = annotation_colors,
    display_numbers = significance_matrix,
    fontsize_number = 13,
    number_color = "blue",
    gaps_col = gaps,
    angle_col = "315",
    na_col = "white",
    background = "transparent",
    breaks = breaks,
    
    # 在这里使用处理过的标签
    labels_row = clean_row_labels, 
    
    # 关键点：设为 NA，这样 pheatmap 会返回对象而不是直接存文件
    filename = NA, 
    
    width = 8.4, height = 7.6,
    legend_position = "bottomright",
    legend = TRUE,
    fontsize_row = 11,
    fontsize_col = 11
  )
  
  # --- 改动 3: 修改绘图对象，将行名设为斜体 ---
  # 找到热图对象中对应 "row_names" (行名) 的部分
  row_label_idx <- which(heatmap$gtable$layout$name == "row_names")
  
  # 将该部分的字体样式 (fontface) 强制改为斜体 ("italic")
  heatmap$gtable$grobs[[row_label_idx]]$gp$fontface <- "italic"
  
  # --- 改动 4: 手动保存修改后的热图 ---
  pdf(save_file, width = 8.4, height = 7.6)
  grid::grid.draw(heatmap$gtable)
  dev.off()
}













############HDP 16S ITS##################################



data_type = "ITS"
transformer = "c"
wt = "wt1"
state = "all"



file_name = paste0('../Anti_data/preprocessed/HDP_ALL_', data_type, '_', transformer, '_', wt, '_AbundInfo0.25_Info_no_antibiotic_as_cov_class3.xlsx')

# HDP_T1 = read.xlsx(file_name, sheet = 1)
# HDP_T2 = read.xlsx(file_name, sheet = 2)
# HDP_T3 = read.xlsx(file_name, sheet = 3)
# HDP_Lmer = read.xlsx(file_name, sheet = 10)
PE_T1 = read.xlsx(file_name, sheet = 4)
PE_T2 = read.xlsx(file_name, sheet = 5)
PE_T3 = read.xlsx(file_name, sheet = 6)
# PE_Lmer = read.xlsx(file_name, sheet = 11)
GH_T1 = read.xlsx(file_name, sheet = 7)
GH_T2 = read.xlsx(file_name, sheet = 8)
GH_T3 = read.xlsx(file_name, sheet = 9)

PEvsGH_T1 = read.xlsx(file_name, sheet = 10)
PEvsGH_T2 = read.xlsx(file_name, sheet = 11)
PEvsGH_T3 = read.xlsx(file_name, sheet = 12)

PE_Pooled = read.xlsx(file_name, sheet = 13)
GH_Pooled = read.xlsx(file_name, sheet = 14)
PEvsGH_Pooled = read.xlsx(file_name, sheet = 15)

# HDP_T1$source <- "HDP_T1"
# HDP_T2$source <- "HDP_T2"
# HDP_T3$source <- "HDP_T3"
PE_T1$source <- "PE_T1"
PE_T2$source <- "PE_T2"
PE_T3$source <- "PE_T3"
GH_T1$source <- "GH_T1"
GH_T2$source <- "GH_T2"
GH_T3$source <- "GH_T3"
PEvsGH_T1$source <- "PE vs. GH (T1)"
PEvsGH_T2$source <- "PE vs. GH (T2)"
PEvsGH_T3$source <- "PE vs. GH (T3)"

PE_Pooled$source <- "Pooled (PE)"
GH_Pooled$source <- "Pooled (GH)"
PEvsGH_Pooled$source <- "Pooled (PE vs. GH)"




if(state=="T1"){
  combined_data <- bind_rows(HDP_T1, PE_T1, GH_T1)
}
if(state=="T23"){
  combined_data <- bind_rows(HDP_T2, PE_T2, GH_T2, HDP_T3,  PE_T3, GH_T3)
}
if(state=="all"){
  combined_data <- bind_rows(PE_T1, GH_T1, PEvsGH_T1, PE_T2, GH_T2, PEvsGH_T2, PE_T3, GH_T3, PEvsGH_T3, PE_Pooled, GH_Pooled, PEvsGH_Pooled)
}



combined_data <- combined_data %>%
  mutate(data_source = 'ITS')




file_name2 = paste0('../Anti_data/preprocessed/HDP_ALL_Taxonomy_g_c_wt2_AbundInfo0.25_Info_no_antibiotic_as_cov_class3.xlsx')


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

PEvsGH_T1_2 = read.xlsx(file_name2, sheet = 10)
PEvsGH_T1_2$source <- "PE vs. GH (T1)"
PEvsGH_T2_2 = read.xlsx(file_name2, sheet = 11)
PEvsGH_T2_2$source <- "PE vs. GH (T2)"
PEvsGH_T3_2 = read.xlsx(file_name2, sheet = 12)
PEvsGH_T3_2$source <- "PE vs. GH (T3)"

PE_Pooled_2 = read.xlsx(file_name2, sheet = 13)
GH_Pooled_2 = read.xlsx(file_name2, sheet = 14)
PEvsGH_Pooled_2 = read.xlsx(file_name2, sheet = 15)

PE_Pooled_2$source <- "Pooled (PE)"
GH_Pooled_2$source <- "Pooled (GH)"
PEvsGH_Pooled_2$source <- "Pooled (PE vs. GH)"


combined_data2 <- bind_rows(PE_T1_2, GH_T1_2, PEvsGH_T1_2, PE_T2_2, GH_T2_2, PEvsGH_T2_2, PE_T3_2, GH_T3_2, PEvsGH_T3_2, PE_Pooled_2, GH_Pooled_2, PEvsGH_Pooled_2)
combined_data2 <- combined_data2[combined_data2$meta_name != "g93", ]
combined_data2 <- combined_data2 %>%
  mutate(data_source = 'MGG')



# 将第一个文件和第二个文件的合并结果再次合并
combined_data <- bind_rows(combined_data2, combined_data)


combined_data <- combined_data %>%
  mutate(name = ifelse(name == "g__uncultured", meta_name, name))
combined_data <- combined_data %>%
  mutate(beta = pmin(pmax(beta, -3), 3))

combined_data <- combined_data[!grepl('g[0-9]', combined_data$name), ]


g_data <- combined_data %>% filter(data_source == "ITS")
s_data <- combined_data %>% filter(data_source == "MGG")


# 找到 data_source 为 '16S' 且在 HDP_T1, PE_T1, GH_T1 中任意一列 p 值小于 0.05 的 name 顺序
g_priority_names <- g_data %>%
  group_by(name) %>%
  filter(any(source %in% c("PE_T2", "GH_T2", "PEvsGH_T2") & pv_adj_FDR < 0.1)) %>%
  ungroup() %>%
  distinct(name) %>%
  arrange(name) %>%  # 对 name 重新排序
  pull(name)

g_priority_df <- g_data %>%
  filter(name %in% g_priority_names) %>%
  distinct(name, .keep_all = TRUE) %>%  
  select(name, meta_name)


write_xlsx(g_priority_df, "Anti_results/ITS_sig.xlsx")

# 找到剩余的 data_source 为 '16S' 的 name 顺序
g_remaining_names <- g_data %>%
  filter(!(name %in% g_priority_names)) %>%
  distinct(name) %>%
  arrange(name) %>%  # 对 name 重新排序
  pull(name)

# 找到 data_source 为 'ITS' 且在 HDP_T1, PE_T1, GH_T1 中任意一列 p 值小于 0.05 的 name 顺序
s_priority_names <- s_data %>%
  group_by(name) %>%
  filter(any(source %in% c("PE_T1", "GH_T1", "PEvsGH_T1") & pv_adj_FDR < 0.1)) %>%
  ungroup() %>%
  distinct(name) %>%
  arrange(name) %>%  # 对 name 重新排序
  pull(name)

s_priority_df <- s_data %>%
  filter(name %in% s_priority_names) %>%
  distinct(name, .keep_all = TRUE) %>%  
  select(name, meta_name)




# 找到剩余的 data_source 为 'ITS' 的 name 顺序
s_remaining_names <- s_data %>%
  filter(!(name %in% s_priority_names)) %>%
  distinct(name) %>%
  arrange(name) %>%  # 对 name 重新排序
  pull(name)

# 最终的 name 顺序
final_row_names <- c(s_priority_names, s_remaining_names, g_priority_names, g_remaining_names)



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

if(data_type=="ITS"){

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
  unique_phylum <- c("Bacteroidetes", "Firmicutes", "Proteobacteria", 
                     "Actinobacteria", "Ascomycota", "Basidiomycota") ##"Ascomycota", "Basidiomycota", "Mucoromycota"
  phylum_colors <- setNames(colorRampPalette(brewer.pal(6, "Set3"))(length(unique_phylum)), unique_phylum)
}
}

if(data_type=="16S"){
  
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

# final_row_order <- heatmap_result$tree_row$order
# final_row_names <- rownames(beta_data_orig)[final_row_order]
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
if(state=="T1"){
  save_file = paste0('Anti_graph/HDP_', data_type, '_', transformer, '_', wt, '_clip_heatmap_T1.pdf')
}
if(state=="T23"){
save_file = paste0('Anti_graph/HDP_', data_type, '_', transformer, '_', wt, '_clip_heatmap_T23.pdf')
}
if(state=="all"){
  save_file = paste0('Anti_graph/HDP_', data_type, '_', transformer, '_', wt, '_clip_heatmap_as_cov_class3.pdf')
}



if(data_type=='ITS'){
  
  beta_data_orig <- apply(beta_data_orig, c(1, 2), function(x) pmax(pmin(x, 1), -1))

  breaks <- seq(-1, 1, length.out = 101)

  # pheatmap(beta_data_orig, 
  #          cluster_rows = FALSE, 
  #          cluster_cols = FALSE, 
  #          annotation_row = row_annotation, 
  #          annotation_colors = annotation_colors, 
  #          display_numbers = significance_matrix,  # 添加显著性标注
  #          number_color = "blue",
  #          gaps_col = gaps,
  #          angle_col = "315",
  #          na_col = "white",
  #          breaks = breaks,
  #          filename = save_file,
  #          # width = 6.3, height =3.5, # for T1
  #          width = 7.83, height =9.754, # for T23
  #          legend_position = "bottomright",
  #          legend = TRUE)
  
  clean_row_labels <- gsub("^g__", "", rownames(beta_data_orig))
  
  # --- 改动 2: 生成热图对象，但先不保存 (filename = NA) ---
  heatmap <- pheatmap(
    beta_data_orig,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = row_annotation,
    annotation_colors = annotation_colors,
    display_numbers = significance_matrix,
    fontsize_number = 13,
    number_color = "blue",
    gaps_col = gaps,
    angle_col = "315",
    na_col = "white",
    background = "transparent",
    breaks = breaks,
    
    # 在这里使用处理过的标签
    labels_row = clean_row_labels, 
    
    # 关键点：设为 NA，这样 pheatmap 会返回对象而不是直接存文件
    filename = NA, 
    
    width = 7.3, height = 7.5,
    legend_position = "bottomright",
    legend = TRUE,
    fontsize_row = 11,
    fontsize_col = 11
  )
  
  # --- 改动 3: 修改绘图对象，将行名设为斜体 ---
  # 找到热图对象中对应 "row_names" (行名) 的部分
  row_label_idx <- which(heatmap$gtable$layout$name == "row_names")
  
  # 将该部分的字体样式 (fontface) 强制改为斜体 ("italic")
  heatmap$gtable$grobs[[row_label_idx]]$gp$fontface <- "italic"
  
  # --- 改动 4: 手动保存修改后的热图 ---
  pdf(save_file, width = 7.3, height = 7.5)
  grid::grid.draw(heatmap$gtable)
  dev.off()

}

