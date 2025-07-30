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

set.seed(555)
rm(list = ls())

df_all1 = read.xlsx("data/Random_forest_HDP.xlsx") # switch to Random_forest_HDP/PE/GH

ggplot(df_all1, aes(x = FPR, y = TPR, color = Label)) +
  geom_line(size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(title = "ROC curves",
       x = "False positive rate",
       y = "True positive rate") +
  theme_minimal(base_size = 15) +
  theme(legend.position = c(0.7, 0.15), 
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = NA),
        legend.key = element_rect(fill = "white", color = NA))

df_all1 = df_all1 %>% arrange(Label, FPR, TPR)
ggplot(df_all1, aes(x = FPR, y = TPR, color = Label)) +
  geom_line(size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(x = "1-Specificity",
       y = "Sensitivity") +
  theme_minimal(base_size = 15) +
  theme(axis.title = element_text(size = 20), 
        legend.position = c(0.65, 0.15), 
        legend.title = element_blank(),
        axis.text = element_text(size = 20, color = "black"),
        legend.background = element_rect(fill = "white", color = NA),
        legend.key = element_rect(fill = "white", color = NA), 
        legend.text = element_text(size = 13),
        panel.border = element_rect(color = "black", fill = NA, size = 1))












