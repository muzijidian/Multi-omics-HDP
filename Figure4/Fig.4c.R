library(dplyr)
library(ggplot2)
library(ggpubr) 
library(ggsignif)
source('./00_constant_params.R')
source('./myUtils.R')

WestLake_HBC_ITS16S_Abundance <-  function(data, gen, cov_names, method='arcsinZ', prevalence=0.25, micro_type='ITS') {
    data = data[!duplicated(data %>% select(id, period)),] 
    explan_var = data %>% dplyr::select(all_of(cov_names))

    min_abund=1e-4
    gen_rule <- paste('^', gen, '\\d+', sep = "")
    count_tab_gen = data %>% dplyr::select(matches(gen_rule))
    all_micro <- names(count_tab_gen)

    count_tab = decostand(count_tab_gen, 'total', 1)
    micro_sel = select_micro(count_tab_gen, min.abundance=min_abund, min.prevalence=prevalence)
    clr_tab = count_tab %>% select(all_of(micro_sel))
    print(paste0('Abund:', dim(clr_tab)))

    data = cbind(explan_var, clr_tab)
    data <- data %>% 
        mutate_at(micro_sel, ~ . + 1e-10) %>%
        mutate_at(micro_sel, sqrt)
    data <- data %>% mutate_at(micro_sel, asin)
    data <- data %>% mutate_at(micro_sel, z_trans)
    data <- data.frame(data)

    return(list(data=data, micro_sel=micro_sel))
}


THSBC_dataprepare <- function(micro, sig, gen='g', wt='wt3'){
    file <- paste0('../data/THSBC_', micro, '.csv')
    data <- readData(file)
    data$PE <- data$preeclampsia
    data$GH <- data$HDP - data$preeclampsia
    data$prob <- data[[wt]]
    data <-  data  %>%  relocate(PE, GH, .after=HDP)
    cat(table(data$period), '\n')

    if (micro != 'Met'){
        sel_names <- c('id', 'period', 'prob', 'GH', 'PE', 'HDP')
        Abundance <- ITS16S_Abundance(data, gen, sel_names, 'arcsinZ', prevalence=0.25, micro_type=micro) 
        data <- Abundance$data
    }
    data <- filter_antibiotic(data)
    data <- data  %>%  filter(hpt_prep == 0)
    cat('Data no antibiotic used: ', dim(data), '\n') 

    Vn <- ifelse(micro == 'ITS', 'V2', 'V1')
    data <- filter(data, period==Vn)  
    data <- data[!is.na(data$HDP), ]

    sel_names <- c('id', 'Type', sig)
    data_HDP <- filter(data, HDP == 1)  
    data_HDP$Type <- 'HDP'
    data_HDP <- data_HDP  %>%  select(all_of(sel_names))
    
    data <- data  %>%  mutate(
        Type=case_when(
            HDP == 0 ~ 'Non-HDP',
            PE  == 1 ~ 'PE',
            GH  == 1 ~ 'GH'
        )
    )  %>%  select(all_of(sel_names))

    data_THSBC <- rbind(data_HDP, data)
    names(data_THSBC) <- c('DATABASE', 'Type', 'Bac')
    data_THSBC$DATABASE <- '1.THSBC'

    return(data_THSBC)
}


WeBirth_dataprepare <- function(micro, sig, gen='g'){
    file <- paste0('../data/WeBirth_', micro, '.csv')
    data <- readData(file)
    data <- data  %>%  filter(!is.na(HDP))
    data$GH <- data$HDP - data$PE
    data <-  data  %>%  relocate(GH, .after=PE)
    cat(table(data$period), '\n')

    if (micro == 'Met'){
        n_names <- length(names(data))
        micro_name <- names(data)[21:n_names]   
        data <- getDataTrans(data, micro_name)
    } else {
        sel_names <- c('id', 'period')  # 'GH', 'PE', 'HDP'
        Abundance <- WestLake_HBC_ITS16S_Abundance(data, gen, sel_names, 'arcsinZ', prevalence=0.25, micro_type=micro) 
        data <- Abundance$data
    }

    data <- filter_hpt_WeBirth(data)
    data <- data  %>%  filter(pre_hpt == 0)
    data <- filter(data, period=='V1')  

    sel_names <- c('id', 'Type', sig)
    data_HDP <- filter(data, HDP == 1)  
    data_HDP$Type <- 'HDP'
    data_HDP <- data_HDP  %>%  select(all_of(sel_names))
     
    data <- data  %>%  mutate(
        Type=case_when(
            HDP == 0 ~ 'Non-HDP',
            PE  == 1 ~ 'PE',
            GH  == 1 ~ 'GH'
        )
    )  %>%  select(all_of(sel_names))

    data_WeBirth <- rbind(data_HDP, data)
    names(data_WeBirth) <- c('DATABASE', 'Type', 'Bac')
    data_WeBirth$DATABASE <- '2.WeBirth'

    return(data_WeBirth)
}


HBC_dataprepare <- function(micro, sig, gen='g'){
    file <- paste0('../data/HBC_', micro, '_1042.csv')
    data <- readData(file)
    data <- data  %>%  filter(!is.na(HDP))
    data$GH <- data$HDP - data$PE
    data <-  data  %>%  relocate(GH, .after=PE)
    cat(table(data$period), '\n')

    if (micro == 'Met'){
        n_names <- length(names(data))
        micro_name <- names(data)[3:966]   
        data <- getDataTrans(data, micro_name)
    } else {
        sel_names <- c('id', 'period')  # , 'GH', 'PE', 'HDP'
        Abundance <- WestLake_HBC_ITS16S_Abundance(data, gen, sel_names, 'arcsinZ', prevalence=0.25, micro_type=micro) 
        data <- Abundance$data
    }
    
    data <- filter_hpt_HBC(data)
    data <- data  %>%  filter(pre_hpt == 0)
    data <- filter(data, period=='V1')  

    sel_names <- c('id', 'Type', sig)
    data_HDP <- filter(data, HDP == 1)  
    data_HDP$Type <- 'HDP'
    data_HDP <- data_HDP  %>%  select(all_of(sel_names))
     
    data <- data  %>%  mutate(
        Type=case_when(
            HDP == 0 ~ 'Non-HDP',
            PE  == 1 ~ 'PE',
            GH  == 1 ~ 'GH'
        )
    )  %>%  select(all_of(sel_names))

    data_HDP <- rbind(data_HDP, data)
    names(data_HDP) <- c('DATABASE', 'Type', 'Bac')
    data_HDP$DATABASE <- '3.HBC'

    return(data_HDP)
}


IQR_filtered_data <- function(pdata){
    Q1 = quantile(pdata[['Bac']], 0.25)
    Q3 = quantile(pdata[['Bac']], 0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR

    filtered_data = pdata[(pdata[['Bac']] >= lower_bound) & (pdata[['Bac']] <= upper_bound),]

    return(filtered_data)
}


main <- function(micro, codes, micro_name, gen='g', outFile="./vioplot.pdf"){
    data_THSBC = THSBC_dataprepare(micro, codes[1], gen)
    data_WeBirth = WeBirth_dataprepare(micro, codes[2], gen)
    data_HBC = HBC_dataprepare(micro, codes[3], gen)
    data = rbind(data_THSBC, data_WeBirth, data_HBC)

    cat('THSBC: ',dim(data_THSBC), '\n')
    cat('WeBirth: ',dim(data_WeBirth), '\n')
    cat('HBC: ', dim(data_HBC), '\n')

    data_THSBC <- IQR_filtered_data(data_THSBC)
    data_WeBirth <- IQR_filtered_data(data_WeBirth)
    data_HBC <- IQR_filtered_data(data_HBC)
    data_filtered = rbind(data_THSBC, data_WeBirth, data_HBC)

    my_comparisons=list(c('Non-HDP', 'HDP'), c('Non-HDP', 'PE'), c('Non-HDP', 'GH'))

    t_test_result <- c()
    for (database in c('1.THSBC', '2.WeBirth', '3.HBC')){
        t_test_compare <- c()
        DATA <- data  %>%  filter(DATABASE == database)
        print(table(DATA$Type))
        for (compare in my_comparisons){
            t_test <- wilcox.test(DATA$Bac[DATA$Type == compare[1]], DATA$Bac[DATA$Type == compare[2]])
            t_test_compare <- rbind(t_test_compare, t_test$p.value)
        }
        t_test_result <- cbind(t_test_result, t_test_compare)
    }
    colnames(t_test_result) <- c('THSBC', 'WeBirth', 'HBC')
    rownames(t_test_result) <- c('HDP', 'PE', 'GH')
    t_test_result <- data.frame(t_test_result)
    print(t_test_result)

    p1=ggplot(data_filtered,aes(x=Type,y=Bac,fill=Type), outlier.shape=NA)+
        guides(fill=guide_legend(), outlier.shape=NA)+
        labs(y = micro_name)+
        geom_violin(trim = F) + 
        # coord_cartesian(ylim = c(-4, 8)) +
        geom_boxplot(width=0.2,position=position_dodge(0.9), outlier.shape=NA) + 
        stat_compare_means(data=data, comparisons = my_comparisons, label = "p.format", label.y = c(2,2.5,3), method='wilcox.test') +     # 
        # stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method='t.test') +
        scale_fill_manual(values = c("GH" = "#f8c672", "HDP" = "#efb0b3", "Non-HDP" = "#7a94ba", "PE" = "#B5665D")) +
        facet_wrap(~DATABASE, nrow = 1) + 
        theme_bw()+
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=14, face = "plain", colour = 'black'),
              axis.text.y = element_text(angle = 0, hjust = 1, size=14, face = "bold", colour = 'black'),
              axis.title.x = element_blank(),
              axis.title.y = element_text(face = "bold", size = 14, angle = 90),
              panel.grid = element_blank(),
              panel.background = element_blank(), 
              plot.background = element_blank(),
              strip.text = element_text(size = 14, face = "bold", color = "black"), 
              strip.background = element_rect(fill = 'lightblue', color = "black", size = 1) 
        ) +
        scale_x_discrete(limits = c("Non-HDP", "HDP", "PE", "GH"))

    print(p1)
    # ggsave(outFile, p1, width=11, height=5, bg = "transparent")
    print(outFile)
}


codes <- c('s50', 's17', 's50')
micro_name <- 's__Clostridium_fessum'
micro <- 'Tax'
main(micro, codes, micro_name, 's', outFile=paste0("./Fig4c/Fig.4c.pdf"))


codes <- c('g130', 'g123', 'g130')
micro_name <- 'g__Cladosporium'
micro <- 'ITS'
main(micro, codes, micro_name, outFile=paste0("./Fig4c/Fig.4c.pdf"))
