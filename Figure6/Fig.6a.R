
library(ggplot2)
library(readxl)

X5 <- read_excel("../data/X5.125718456_G_PE_withouGH.xlsx")
# X11 <- read_excel("../data/X11.125123745_T_PE_withouGH.xlsx")
# X22 <- read_excel("../data/X22.44647223_T_PE_withouGH.xlsx")
# X4 <- read_excel("../data/X4.99614841_T_PE_withouGH.xlsx")
# By replacing different data, you can get their own graphics
X5$significance <- ifelse(X5$adonis.p == 0.001, "***", 
                          ifelse(X5$adonis.p < 0.01, "**", 
                                 ifelse(X5$adonis.p < 0.05, "*", "")))
p1<-ggplot(X5, aes(y = group, x = adonis.r2, fill = micro_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.3), width = 0.7) + 
  geom_text(aes(label = significance), 
            position = position_dodge(width = 0.5), 
            hjust = -0.2, size = 5, color = "black") +  
  scale_fill_manual(values = c("16S" = "#95C6C6", "Taxonomy_g" = "#7EABCA", "Taxonomy_s" = "#FFB2B9", "ITS" = "#ebdcb2", "Metabolites" = "#472D7B3a")) + 
  labs(title = " ", 
       x = NULL, 
       y = NULL,
       fill = " ") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    text = element_text(size = 18),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black",size = 14), 
    axis.title.x.top = element_text(size = 14, color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x.bottom = element_blank(),
    axis.line.x.top = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.position = "none",  
    legend.justification = c(1, 0), 
    legend.title = element_blank(),  
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = margin(t = 5, r = 25, b = 5, l = 5)
  ) +
  scale_x_continuous(position = "top", 
                     limits = c(0, 0.004),  
                     expand = c(0, 0)) + 
  facet_wrap(~ micro_type, scales = "free_y", ncol = 1)

ggsave(
  filename = "X5_PE_cor_new.tiff", 
  plot = p1,             
  width = 6,                     
  height = 3,                   
  units = "in",                  
  dpi = 1000,                    
  device = "tiff"               
)
