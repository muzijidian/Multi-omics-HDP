library(pheatmap)
library(openxlsx)
library(tidyr)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(stringr)
library(readxl)
library(Boruta)
library(RColorBrewer)
library(viridis)
library(colorspace)
library(randomForest)
library(caret)
library(pROC)

setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code')
source('../R_code/myUtils.R')
set.seed(555)
rm(list = ls())



plot.importance.final = function(data){
  name <- deparse(substitute(data))


  data$Variable <- gsub("`|\\d+$", "", data$Variable)
  data$Variable <- gsub("^meta_", "", data$Variable)
  data$Variable <- gsub(" Acid", " acid", data$Variable)


  ggplot(data, aes(x = reorder(Variable, Importance), y = Importance, fill = Group)) +
    geom_bar(stat = "identity", width = 0.7) + 
    coord_flip() + 
    labs(title = "Variable importance in random forest",
         x = " ",
         y = "Importance") +
    theme_minimal(base_size = 12) + 
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
      axis.text = element_text(size = 9, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      panel.grid.major = element_line(color = "grey80"),  
      panel.grid.minor = element_blank(),  
      panel.background = element_blank(), 
      legend.position = "none"
    ) +
    scale_fill_manual(values = c(
      "16S" = "#95C6C6",
      "MetaG_g" = "#7EABCA",
      "MetaG_s" = "#FFB2B9",
      "ITS" = "#ebdcb2",
      "Metabolites" = "#472D7B3a",
      "TFS" = "steelblue"
    ))  


  filename = paste0('./graph/', name, '_wt_batch_HDP.tiff')
  ggsave(filename,
         width = 7,
         height = 9,
         dpi = 800,
         compression = "lzw")
}

fina_data_all = read.xlsx("../data/roc_o_importance_HDP.xlsx")
plot.importance.final(fina_data_all)


