# setwd('/data/TalenWang/Codes/myDepression/src/')
library(dplyr)
library(openxlsx)
library(reshape2)
library(stringr)
library(foreign)
library(haven)
library(vegan)


readData <- function(file) {
    file_split <- strsplit(file, '\\.')[[1]]
    file_suffix <- file_split[length(file_split)] 
    print(paste0('Reading file: ', file))
    if (file_suffix == 'xlsx'){
        data <- read.xlsx(file)
    }
    else if (file_suffix == 'dta') {
       data <- read_dta(file)
    } 
    else if (file_suffix == 'csv') {
        data <- read.csv(file)
    }
    else if (file_suffix == 'rds') {
        data <- readRDS(file)
    }
    else {
        print('Dont support such file!!!')
        q()
    }

    return(data)
}

DataCovPreprocess <- function(data){
    data <- data %>% mutate(
        smk = case_match(smk, 1:4~1, 5~0, .default = 0),
        drk = case_match(drk, 1~0, 2:3~1, .default = 0),
        parity = case_match(parity, 0~0, 1:2~1, .default = 0),
        edu = case_match(edu, 1:4~0, 5:8~1, 9~0, .default = 0),
    )
    return(data)
}

getDataBMIandwk <- function(data) {
    data = data %>% mutate(wk_b = wk_b_wk_pad + wk_b_day_pad/7,
                       wk_c = wk_c_wk_pad + wk_c_day_pad/7,
                       wk_delivery = wk_delivery_wk + wk_delivery_day/7,
                       BMI_prep = wt_prep/height_a/height_a*10000,
                       BMI_a = weight_a/height_a/height_a*10000,
                       BMI_b = weight_b/height_a/height_a*10000,
                       BMI_c = weight_c/height_a/height_a*10000)
    return(data)
}


getDataTrans <- function(data, meta_names){
    data <- data %>% mutate_at(meta_names, log) %>% 
        mutate_at(meta_names, scale)

    return(data)
}


getDataEPDSbool <- function(data, thresh=13){
    data <- data %>% mutate(
        epds_score_bool = case_match(epds_score, -10~NA, 0:thresh-1~0, thresh:30~1)) 
    return(data)
}


DataMelt <- function(datac){
    tab_wk <- datac %>% select(id, period, wk_a, wk_b, wk_c) %>%
                mutate(wk = case_when(
                                period == "V1" ~ wk_a,
                                period == "V2" ~ wk_b,
                                period == "V3" ~ wk_c
                        )) %>% select(id, period, wk) 
    print(paste0('tab_wk: ', dim(tab_wk)))
    # print(head(tab_wk))

    tab_BMI <- datac %>% select(id, period, BMI_a, BMI_b, BMI_c) %>%
                mutate(BMI = case_when(
                                period == "V1" ~ BMI_a,
                                period == "V2" ~ BMI_b,
                                period == "V3" ~ BMI_c
                        )) %>% select(id, period, BMI) 
    print(paste0('tab_BMI: ', dim(tab_BMI)))
    # print(head(tab_BMI))

    tab_weight <- datac %>% select(id, period, weight_a, weight_b, weight_c) %>%
                mutate(weight = case_when(
                                period == "V1" ~ weight_a,
                                period == "V2" ~ weight_b,
                                period == "V3" ~ weight_c
                        )) %>% select(id, period, weight) 
    print(paste0('tab_weight: ', dim(tab_weight)))
    # print(head(tab_weight))

    tab_epds <- datac %>% select(id, period, epds_a_total, epds_b_total, epds_c_total) %>%
            mutate(epds_score = case_when(
                            period == "V1" ~ epds_a_total,
                            period == "V2" ~ epds_b_total,
                            period == "V3" ~ epds_c_total
                    )) %>% select(id, period, epds_score)
    tab_epds <- getDataEPDSbool(tab_epds)
    tab_epds$epds_score[tab_epds$epds_score == -10] <- NA
    print(paste0('tab_epds: ', dim(tab_epds)))

    # 必须要加这一列，才能对齐（不知道为啥）
    tab_wk$X <- c(1:length(tab_wk[, 'id']))
    tab_BMI$X <- c(1:length(tab_BMI[, 'id']))
    tab_weight$X <- c(1:length(tab_weight[, 'id']))
    tab_epds$X <- c(1:length(tab_epds[, 'id']))
    # tab_meta$X <- c(1:length(tab_meta[, 'id']))

    tab_long <- tab_wk %>% left_join(tab_BMI) %>% left_join(tab_weight) %>% left_join(tab_epds) %>% select(-X)
    # %>% left_join(tab_meta)
    # print(paste0('tab_long: ', dim(tab_long)))

    print(paste0('tab_long: ', dim(tab_long)))
    print(names(tab_long))

    # write.csv(tab_long, save_data_path)

    return(tab_long)
}


get_Base_cov <- function() {
    prefile <- '../../phenotype_tsbc_7000.dta'
    pre_data <- readData(prefile)
    cov_names <- c('id', 'epds_a_total', 'epds_b_total', 'epds_c_total', 'age', 'smk', 'drk', 'parity', 'edu', 'height_a', 'weight_a', 'weight_b', 'weight_c', 'wt_prep', 'wk_a', 'wk_b_wk_pad', 'wk_b_day_pad', 'wk_c_wk_pad', 'wk_c_day_pad', 'wk_delivery_wk', 'wk_delivery_day', 'sbp_a_1','sbp_a_2','sbp_b_1','sbp_b_2','sbp_c_1','sbp_c_2','sbp_ave','dbp_a_1','dbp_a_2','dbp_b_1','dbp_b_2','dbp_c_1','dbp_c_2','dbp_ave')
    cov_data <- pre_data %>% select(all_of(cov_names))
    print(paste('cov_data: ', dim(cov_data)))

    cov_data <- getDataBMIandwk(cov_data)
    cov_data <- DataCovPreprocess(cov_data)
    cov_data <- cov_data %>% mutate_at(c('age', 'BMI_prep'), Hmisc::impute, fun=mean) %>% 
    mutate_at(c('smk', 'drk', 'edu'), Hmisc::impute, 0)

    write.csv(cov_data, '../data/Base_cov_info.csv', row.names=F)
} 


filter_antibiotic <- function(data) {
    data_no_antibiotic <- readData('../data/THSBC_merged_data.xlsx')
    data_no_antibiotic$id <- as.numeric(data_no_antibiotic$id)
    data_no_antibiotic <- data_no_antibiotic  %>%  select(c('id', 'diabetes_prep','gdm_ever','ghpt_ever','family_diabetes','family_hpt', 'aspirin_painkiller_use', 'no_antibiotic_use_a', 'no_antibiotic_use_b', 'no_antibiotic_use_c', 'hpt_prep'))

    data <- data  %>%  left_join(data_no_antibiotic, by='id')
    data <- data  %>%  mutate(
        no_antibiotic=case_when(period=='V1' ~ no_antibiotic_use_a,
                                period=='V2' ~ no_antibiotic_use_b,
                                period=='V3' ~ no_antibiotic_use_c,
        ),
        no_antibiotic_use_any=ifelse((no_antibiotic_use_a+no_antibiotic_use_b+no_antibiotic_use_c)>0, 1, 0))
    # data <- data  %>%  filter(no_antibiotic == 1)
    # data <- data  %>%  filter(no_antibiotic_use_any == 1)

    return(data)
}


filter_hpt_HBC <- function(data) {
    data_new_covs <- readData('../data/HBC_final_covariates.csv')
    data_new_covs$id <- as.numeric(data_new_covs$id)
    data_new_covs <- data_new_covs  %>%  select(c('id', 'age','BMI_prep','pa_cat','smk_cat','drk_cat','edu_cat', 'hpt_family','hdp_ever','gdm_ever','dia_family', 'pre_hpt', 'HDP', 'PE'))
    data_new_covs$GH <- data_new_covs$HDP - data_new_covs$PE

    data <- data  %>%  left_join(data_new_covs, by='id')
    # data <- data  %>%  filter(ab_type_0 == 1)

    return(data)
}

filter_hpt_WeBirth <- function(data) {
    data_new_covs <- readData('../data/WeBirth_merged_data.xlsx')
    data_new_covs$id <- as.character(data_new_covs$id)
    data_new_covs <- data_new_covs  %>%  select(c('id', "antibiotic_use", "pre_dia", "pre_hpt", "pre_hdp", "pre_gdm", "dia_family", "hpt_family",'HDP','PE'))
    data_new_covs$GH <- data_new_covs$HDP - data_new_covs$PE

    data <- data  %>%  left_join(data_new_covs, by='id')
    # data <- data  %>%  filter(ab_type_0 == 1)

    return(data)
}


select_micro <- function(x, min.abundance=1e-4, min.prevalence=0.05){
  # 输入为count_tab或relab_tab
  abund_mean = colMeans(decostand(x, 'total', 1))
  prevalence = colMeans(x > 0)
  micro_sel = names(which(abund_mean > min.abundance & prevalence > min.prevalence))

  return(micro_sel)
}


ITS16S_Abundance <-  function(data, gen, cov_names, method='arcsinZ', prevalence=0.25, micro_type='ITS') {
    # data <- filter(data, !is.na(prob))
    # 去重
    # data = data[!duplicated(data %>% select(id, period)),] 
    # 读取EPDS诊断和部分协变量的数据
    explan_var = data %>% dplyr::select(all_of(cov_names))

    # 选择genus水平数据，g开头的列
    min_abund=1e-5
    if (gen == 'path'){min_abund=1e-5}

    gen_rule <- paste('^', gen, '\\d+', sep = "")
    count_tab_gen = data %>% dplyr::select(matches(gen_rule))
    all_micro <- names(count_tab_gen)

    if (method == c('clr')){
        micro_sel = select_micro(count_tab_gen, min.prevalence=prevalence) # 极低丰度high
        count_tab = count_tab_gen %>% select(all_of(micro_sel))
        clr_tab = decostand(count_tab, 'clr', 1, pseudocount = 1) # clr
    } else {
        count_tab = decostand(count_tab_gen, 'total', 1)
        micro_sel = select_micro(count_tab_gen, min.abundance=min_abund, min.prevalence=prevalence) # 极低丰度high
        clr_tab = count_tab %>% select(all_of(micro_sel))

        # print('cccccccccc')
    } 

    # 合并数据
    data = cbind(explan_var, clr_tab)
    # data_names <- names(data)
    # print(head(data, n=1))

    # print(paste0('CHeck:', max(data$g50)))

    # arcsinZ和sqrtZ的后续变换
    if (method == c('arcsinZ')) {
        data <- data %>% 
          mutate_at(micro_sel, ~ . + 1e-10) %>%
          mutate_at(micro_sel, sqrt)

        data <- data %>% mutate_at(micro_sel, asin)
        data <- data %>% mutate_at(micro_sel, z_trans)
        # data <- data %>% mutate_at(micro_sel, ~ ifelse(. > 4 | . < -4, NA, .))

        # print('zzzzzzzzz')

    } else if (method == c('sqrtZ')) {
        data <- data %>% 
          mutate_at(micro_sel, ~ . + 1e-10) %>%
          mutate_at(micro_sel, sqrt)
        data <- data %>% mutate_at(micro_sel, z_trans)
        # data <- data %>% mutate_at(micro_sel, ~ ifelse(. > 4 | . < -4, NA, .))
    }

    data <- filter(data, !is.na(prob))
    data <- data.frame(data)
    # print(paste0('CHeck:', max(data$g50)))

    # Abund Info
    drop_micro = setdiff(all_micro, micro_sel)
    micro_names = c(micro_sel, drop_micro)

    micro_sel_m = rep(method, length(micro_sel))
    drop_micro_m = rep('drop', length(drop_micro))
    micro_method = c(micro_sel_m, drop_micro_m)

    prevalence_cols = colMeans(count_tab_gen > 0)
    prevalence_cols = prevalence_cols[micro_names]

    Abund_Info <- data.frame(
        micro = micro_names, 
        method = micro_method, 
        prevalence = prevalence_cols)

    Abund_file = paste0('./violin_plot/', micro_type, '_',gen, '_', method, '_prevalence_', prevalence, '.csv')
    write.csv(Abund_Info, file=Abund_file, row.names=F)

    return(list(data=data, micro_sel=micro_sel))
}


# ITS16S_Abundance <-  function(data, gen, cov_names, method='clr', prevalence_low=0.15, prevalence_high=0.20, micro_type='ITS') {
#     # data <- filter(data, !is.na(prob))
#     # 去重
#     data = data[!duplicated(data %>% select(id, period)),] 
#     # 读取EPDS诊断和部分协变量的数据
#     explan_var = data %>% dplyr::select(all_of(cov_names))

#     # 选择genus水平数据，g开头的列
#     gen_rule <- paste('^', gen, '\\d+', sep = "")
#     count_tab_gen = data %>% dplyr::select(matches(gen_rule))
#     all_micro <- names(count_tab_gen)

#     if (method == c('clr')){
#         micro_sel_low = select_micro(count_tab_gen, min.prevalence=prevalence_low) # 极低丰度low
#         micro_sel_high = select_micro(count_tab_gen, min.prevalence=prevalence_high) # 极低丰度high
#         count_tab = count_tab_gen %>% select(all_of(micro_sel_high))
#         clr_tab = decostand(count_tab, 'clr', 1, pseudocount = 1) # clr
#     } else {
#         count_tab = decostand(count_tab_gen, 'total', 1)
#         micro_sel_low = select_micro(count_tab, min.prevalence=prevalence_low) # 极低丰度low
#         micro_sel_high = select_micro(count_tab_gen, min.prevalence=prevalence_high) # 极低丰度high
#         clr_tab = count_tab %>% select(all_of(micro_sel_high))
#     } 

#     # 低丰度菌群处理，所有prevalence~0.1内的有效数值的置1，0保持为0
#     low_abund_micro <- setdiff(micro_sel_low, micro_sel_high)
#     low_abund_data <- count_tab_gen  %>%  select(all_of(low_abund_micro))
#     low_abund_data <- low_abund_data  %>%  mutate(across(everything(), ~ifelse(. > 0, 1, 0)))

#     # Abund Info
#     drop_micro = setdiff(all_micro, micro_sel_low)
#     micro_names = c(micro_sel_high, low_abund_micro, drop_micro)

#     micro_sel_m = rep(method, length(micro_sel_high))
#     low_abund_m = rep('01', length(low_abund_micro))
#     drop_micro_m = rep('drop', length(drop_micro))
#     micro_method = c(micro_sel_m, low_abund_m, drop_micro_m)

#     prevalence = colMeans(count_tab_gen > 0)
#     prevalence = prevalence[micro_names]

#     Abund_Info <- data.frame(
#         micro = micro_names, 
#         method = micro_method, 
#         prevalence = prevalence)

#     # print(head(Abund_Info, n=3))
#     # stop()

#     write.csv(Abund_Info, paste0('../results/', micro_type, '_', method,'_prevalence', prevalence_low, '-', prevalence_high, 'set01.csv'), row.names=F)
#     # stop()

#     # 合并数据
#     data = cbind(explan_var, clr_tab, low_abund_data)
#     data <- filter(data, !is.na(prob))
#     data_names <- names(data)
#     print(head(data, n=1))

#     # arcsinZ和sqrtZ的后续变换
#     if (method == c('arcsinZ')) {
#         data <- data %>% mutate_at(micro_sel_high, sqrt)
#         data <- data %>% mutate_at(micro_sel_high, asin)
#         data <- data %>% mutate_at(micro_sel_high, z_trans)
#         data <- data %>% mutate_at(micro_sel_high, ~ ifelse(. > 4 | . < -4, NA, .))

#     } else if (method == c('sqrtZ')) {
#         data <- data %>% mutate_at(micro_sel_high, sqrt)
#         data <- data %>% mutate_at(micro_sel_high, z_trans)
#         data <- data %>% mutate_at(micro_sel_high, ~ ifelse(. > 4 | . < -4, NA, .))
#     } else {
#         # print('CLR micro filter')
#         select_micro_data <- data  %>%  select(all_of(micro_sel_high))
#         mean_select <- colMeans(select_micro_data)
#         std_select <- sapply(select_micro_data, sd)
#         range_selsct_low <- mean_select - 4*std_select
#         range_selsct_high <- mean_select + 4*std_select

#         for (microi in micro_sel_high){
#             data <- data %>% mutate_at(microi, ~ ifelse(. > range_selsct_high[microi] | . < range_selsct_low[microi], NA, .))
#         }
        
#         # print(colSums(is.na(data)))
#         # stop()

#     }

#     data <- data.frame(data)
#     return(list(data=data, micro_sel=micro_sel_high))
# }

z_trans <- function(x){
    (x-mean(x, na.rm=T))/sd(x, na.rm = T)
}


ITS16S_Abundance_AsinZ <-  function(data, gen, cov_names, method='c', prevalence_low=0.10, prevalence_high=0.25, micro_type='ITS') {
    # data <- filter(data, !is.na(prob))
    # 去重
    data = data[!duplicated(data %>% select(id, period)),] 
    # 读取EPDS诊断和部分协变量的数据
    explan_var = data %>% dplyr::select(all_of(cov_names))

    # 选择genus水平数据，g开头的列
    min_abund=1e-4
    # if (gen == 'path'){min_abund=1e-5}

    gen_rule <- paste('^', gen, '\\d+', sep = "")
    count_tab_gen = data %>% dplyr::select(matches(gen_rule))
    all_micro <- names(count_tab_gen)

    count_tab = decostand(count_tab_gen, 'total', 1)
    micro_sel_low = select_micro(count_tab, min.abundance=min_abund, min.prevalence=prevalence_low) # low
    micro_sel_high = select_micro(count_tab, min.abundance=min_abund, min.prevalence=prevalence_high) # high

    if (method == c('c')){
        # c型菌群处理，只取大于高阈值的菌进行AsinZ变换
        micro_sel = micro_sel_high
        clr_tab = count_tab %>% select(all_of(micro_sel))
      

        # 数据中的0替换为NA
        # clr_tab <- clr_tab  %>% mutate_at(micro_sel, ~ ifelse(. == 0, NA, .))

        # arcsinZ的后续变换
        clr_tab <- clr_tab %>% 
          mutate_at(micro_sel, ~ . + 1e-10) %>%
          mutate_at(micro_sel, sqrt)
        clr_tab <- clr_tab %>% mutate_at(micro_sel, asin)
        clr_tab <- clr_tab %>% mutate_at(micro_sel, z_trans)
        # clr_tab <- clr_tab %>% mutate_at(micro_sel, ~ ifelse(. > 4 | . < -4, NA, .))

    } else {
        # b型菌群处理，所有prevalence_low~high内的有效数值的置1，0保持为0
        micro_sel = setdiff(micro_sel_low, micro_sel_high)
        clr_tab = count_tab %>% select(all_of(micro_sel))

        # 01变换
        clr_tab <- clr_tab  %>% mutate(across(everything(), ~ifelse(. > 0, 1, 0)))
    } 

    # Abund Info
    low_abund_micro <- setdiff(micro_sel_low, micro_sel_high)
    drop_micro = setdiff(all_micro, micro_sel_low)
    micro_names = c(micro_sel_high, low_abund_micro, drop_micro)

    micro_sel_m = rep('AsinZ', length(micro_sel_high))
    low_abund_m = rep('01', length(low_abund_micro))
    drop_micro_m = rep('drop', length(drop_micro))
    micro_method = c(micro_sel_m, low_abund_m, drop_micro_m)

    prevalence = colMeans(count_tab_gen > 0)
    prevalence = prevalence[micro_names]

    Abund_Info <- data.frame(
        micro = micro_names, 
        method = micro_method, 
        prevalence = prevalence)

    write.csv(Abund_Info, paste0('../results/HDP/', micro_type, '_', gen, '_', method,'_prevalence', prevalence_low, '-', prevalence_high, 'set01.csv'), row.names=F)
    # stop()


    # 合并数据
    data = cbind(explan_var, clr_tab)
    data <- filter(data, !is.na(prob))
    data <- data.frame(data)

    # print(head(Abund_Info, n=3))
    # print(head(data, n=3))
    # stop()

    return(list(data=data, micro_sel=micro_sel))
}