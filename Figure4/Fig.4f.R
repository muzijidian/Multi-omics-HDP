library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)


Draw_Inter_BOX <- function(data_folder, inter_micro, micro, micro_name, target, Vn, save_file, wt='wt1'){
  file <- paste0(data_folder, 'Interaction_Logistic_ALL_', inter_micro, '_', micro, '_', target, '_', Vn, '_high_', wt, '.csv')
  data_h <- read.csv(file)
  data_h <- data_h  %>%  filter(inf_flag == 0)

  file <- paste0(data_folder, 'Interaction_Logistic_ALL_', inter_micro, '_', micro, '_', target, '_', Vn, '_low_', wt, '.csv')
  data_l <- read.csv(file)
  data_l <- data_l  %>%  filter(inf_flag == 0)

  file <- data_file[inter_micro][[1]]
  data_o <- read.csv(file)
  gen_rule <- ifelse(inter_micro=='Taxonomy_s', '^s\\d+', '^g\\d+')
  data_relAbund <- data_o  %>% select(matches(gen_rule))
  mean_relAbund <- sapply(data_relAbund, mean)

  plotData <- data_h %>% select('meta_name', 'h_pv'='pv', 'h_beta'='beta')
  plotData <- plotData  %>% inner_join(data_l %>% select('meta_name', 'l_pv'='pv', 'l_beta'='beta'), by='meta_name')

  plotData <- plotData  %>%  filter(meta_name != sub(".*_", "", micro))

  plotData <- plotData %>% mutate(Neis_high= -sign(h_beta)*(log10(h_pv)))
  plotData <- plotData %>% mutate(Neis_low= -sign(l_beta)*(log10(l_pv)))

  plotData$relAbund <- sapply(plotData$meta_name, function(x) mean_relAbund[names(mean_relAbund) == x])

  plotData$relAbund <- plotData$relAbund / max(plotData$relAbund)

  plotData <- plotData  %>%  filter(abs(h_beta)<=5)  %>%  filter(abs(l_beta)<=5)

  for(i in c(1:nrow(plotData))){
    if(abs(plotData$Neis_high[i]) < -log10(0.05)){
      if(abs(plotData$Neis_low[i]) < -log10(0.05)){
        plotData$Significance[i] <- "Non-signif"
        plotData$color[i] <- "gray"
      }else{
        plotData$Significance[i] <- "Sig_in_Low"
        plotData$color[i] <- "indianred"
      }
    }else{
      if(abs(plotData$Neis_low[i]) < -log10(0.05)){
        plotData$Significance[i] <- "Sig_in_High"
        plotData$color[i] <- "skyblue3"
      }else{
        plotData$Significance[i] <- "Sig_in_Both"
        plotData$color[i] <- "gold3"
      }
    }
  }

  Micro_Info <- read.xlsx(Info_file[inter_micro][[1]])
  Micro_Info <- Micro_Info  %>%  select(ID, name)
  plotData <- plotData  %>%  left_join(Micro_Info, by=c('meta_name'='ID'))
  plotData <- plotData  %>%  mutate(name=ifelse(Significance=='Non-signif', NA, name))  

  plotData$name <- ifelse(grepl("^g\\d+$", plotData$name), NA, plotData$name)  
  plotData$name <- ifelse(plotData$name=='g__uncultured', NA, plotData$name)   

  # x_max <- max(abs(plotData$Neis_high))
  # y_max <- max(abs(plotData$Neis_low))
  x_max <- 3
  y_max <- 3

  plotData <- plotData  %>%  arrange(Neis_high, Neis_low)

  # plotting ================================
  mytheme = theme_classic(base_line_size = 1) +
              theme(
                plot.title = element_text(size = 15, color = "black", hjust = 0.5),
                strip.background = element_rect(fill = "#f4e1da", color = 'black'),
                axis.title = element_text(size = 15, color = "black", face = "bold", vjust = 1.9, hjust = 0.5), 
                axis.text.x = element_text(size = 13, color = "black", face = "bold", angle = 20, hjust = 1, vjust = 1), 
                axis.text.y = element_text(size = 15, color = "black", face = "plain"),
                strip.text = element_text(size = 14, face = "bold"), 
                # strip.background = element_rect(colour = "black", fill = NA, size = 1.5),
                panel.border = element_rect(colour = "black", fill = NA, size = 1.5), 
                axis.line = element_blank()
              )

  # View(plotData)
  dx=5

  p <- ggplot(plotData, aes(x=Neis_high, y=Neis_low, label=name)) +
    geom_point(aes(color=Significance, size=relAbund)) +
    # coord_cartesian(xlim = c(-x_max, x_max), ylim = c(-y_max, y_max)) +
    scale_color_manual(values = c("Non-signif"="gray","Sig_in_Low"="indianred","Sig_in_High"="skyblue3","Sig_in_Both"="gold3"))+
    # scale_color_manual(values=plotData$color) +
    scale_size(range = c(1,3)) +
    # geom_text_repel(aes(label=taxon)) +
    theme_bw() + 
    theme(legend.position="none",
          panel.grid = element_blank(), 
          panel.background = element_blank(), 
          plot.background = element_blank(),
          axis.text.x = element_text(size = 12, color = "black", face = "plain"),
          axis.text.y = element_text(size = 12, color = "black", face = "plain"),
          axis.title.x = element_text(face = "plain", size = 12, angle = 0),
          axis.title.y = element_text(face = "plain", size = 12, angle = 90)) +
    geom_vline(xintercept = c(-1.30103,1.30103), linetype="dashed", color="darkgray") +
    geom_hline(yintercept = c(-1.30103,1.30103), linetype="dashed", color="darkgray") +
    labs(x = bquote('Directionality' ~ "×" ~ '-log('~italic("P")~')' ~ .(micro_name) ~ 'high'),  
         y = bquote('Directionality' ~ "×" ~ '-log('~italic("P")~')' ~ .(micro_name) ~ 'low')) +
    geom_text(aes(label=name), size=4, check_overlap = T)  +
    scale_x_continuous(limits = c(-dx, dx), breaks = seq(-dx, dx, by = 1), labels = seq(-dx, dx, by = 1)) +
    scale_y_continuous(limits = c(-dx, dx), breaks = seq(-dx, dx, by = 1), labels = seq(-dx, dx, by = 1))

  ggsave(p, filename=save_file, device = "pdf", width = 6, height = 5)
}


Inter_Plot <- function(target, data_folder, micro_list, micro_name_list, save_folder){
  for (inter_micro in c('16S', 'ITS', 'Taxonomy_g', 'Taxonomy_s')){
    for (micro in micro_list){
      micro_name <- micro_name_list[micro][[1]]
      save_file <- paste0(save_folder, target, '_', inter_micro, '_', micro, '_', micro_name, '_V2_pv.pdf')
      Draw_Inter_BOX(data_folder, inter_micro, micro, micro_name, target, 'V2', save_file)
      print(paste0(save_file, ' has done!'))
    }
  }
}

Info_16S <- '../Info_16S_Batch123_phylum.xlsx'
Info_ITS <- '../Info_ITS_Batch123_phylum.xlsx'
Info_TAX <- '../Info_Taxonomy_phylum.xlsx'
Info_file <- list('16S'=Info_16S, 'ITS'=Info_ITS, 'Taxonomy_g'=Info_TAX, 'Taxonomy_s'=Info_TAX)

Data_16S <- './Weighted_ALL_16S.csv'
Data_ITS <- './Weighted_ALL_ITS.csv'
Data_TAX <- './Weighted_ALL_Taxonomy.csv'
data_file <- list('16S'=Data_16S, 'ITS'=Data_ITS, 'Taxonomy_g'=Data_TAX, 'Taxonomy_s'=Data_TAX)

#####################################################HDP Micros##############################################
target <- 'HDP'
data_folder <- './Fig4f/'

micro_list <- c('Taxonomy_s42', 'Taxonomy_s50', 'Taxonomy_s152')
micro_name_list <- list('Taxonomy_s42'='s__Clostridium_sp_AM49_4BH', 
                        'Taxonomy_s50'='s__Clostridium_fessum', 
                        'Taxonomy_s152'='s__Clostridium_sp_AM22_11AC')

save_folder <- './Fig4f/'
Inter_Plot(target, data_folder, micro_list, micro_name_list, save_folder)
