rm(list = ls())
setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code')
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
library(reshape2)
library(mRMRe)
library(PRROC)
library(patchwork)
library(doParallel)
source('help_func.R')

source('../R_code/myUtils.R')

# ==============================================================================
# 设置全局种子 - 整个run使用同一个种子
# ==============================================================================
GLOBAL_SEED <- 555  # 可以修改这个数字
set.seed(GLOBAL_SEED)

# --- 读取特征列表 ---
f_met <- read_excel('Anti_results/Met_sig.xlsx')[['meta_name']]
f_its <- read_excel('Anti_results/ITS_sig.xlsx', sheet = 1)[['meta_name']]
f_bac <- read_excel('Anti_results/Taxonomy_s_sig.xlsx', sheet = 1)[['meta_name']]


# --- 读取数据 ---
met_dat <- read.csv('data/Weighted_ALL_Met_516.csv', check.names = F)%>% 
  filter(period == 'V1')
its_dat <- read.csv('data/Weighted_ALL_ITS.csv', check.names = F)%>% 
  filter(period == 'V1')
bac_dat <- read.csv('data/Weighted_ALL_Taxonomy.csv', check.names = F)%>% 
  filter(period == 'V1')

met_dat <- filter(met_dat, !is.na(wt2))
met_dat <- met_dat %>% 
  filter(!is.na(HDP)) %>% 
  mutate(
    group = case_when(
      HDP == 1 & preeclampsia == 1 ~ "PE",
      HDP == 1 & preeclampsia == 0 ~ "GH",
      HDP == 0                    ~ "NP",
      TRUE                        ~ NA_character_
    )
  )

met_dat$id <- as.numeric(met_dat$id)
its_dat$id <- as.numeric(its_dat$id)
bac_dat$id <- as.numeric(bac_dat$id)

cor_list = c("id", "period", "age", "parity", "edu", "smk", "drk", "wk", "BMI_prep", "group")
cov_data <- met_dat %>%
  select(all_of(cor_list))

cov_names <- setdiff(names(cov_data), c('id', 'group', 'period'))

data_all <- cov_data %>%
  left_join(met_dat %>% select(id, all_of(f_met)), by = 'id') %>%
  left_join(its_dat %>% select(id, all_of(f_its)), by = 'id') %>%
  left_join(bac_dat %>% select(id, all_of(f_bac)), by = 'id') %>%
  tidyr::drop_na()

data_all <- data_all %>% 
  mutate(
    smk = as.factor(smk),
    drk = as.factor(drk),
    edu = as.factor(edu),
    parity = as.factor(parity),
    group = as.factor(group)
  )

ft_combos <- list(
  Baseline = cov_names,
  Metabolome = f_met,
  Fungi = f_its,
  Bacteria = f_bac,
  `Met+Fun` = c(f_met, f_its),
  `Met+Bac` = c(f_met, f_bac),
  `Fun+Bac` = c(f_its, f_bac),
  `Multi-Omics` = c(f_met, f_its, f_bac)
)

tasks <- list(
  GH_vs_NP = list(filter = c('NP', 'GH'), case = 'GH', ctrl = 'NP'),
  Combined_vs_NP = list(filter = c('NP', 'GH', 'PE'), case = c('GH', 'PE'), ctrl = 'NP'),
  PE_vs_NP = list(filter = c('NP', 'PE'), case = 'PE', ctrl = 'NP')
)

# ==============================================================================
# 移除 manual_seed_map，使用全局种子生成策略
# ==============================================================================

# --- 初始化结果容器 ---
final_auc_summary <- data.frame() 
final_roc_tab <- data.frame()
final_imp_tab <- data.frame()

cl <- makeCluster(detectCores() - 2)
registerDoParallel(cl)

target_fpr_grid <- seq(0, 1, length.out = 1001)
target_spec_grid <- 1 - target_fpr_grid

for (task_name in names(tasks)) {
  
  task <- tasks[[task_name]]
  cat(sprintf('\n#### Processing Task: %s ####\n', task_name))
  
  curr_data <- data_all %>% filter(group %in% task$filter)
  curr_data$y <- factor(ifelse(curr_data$group %in% task$case, 'x1', 'x0'), levels = c('x0', 'x1'))
  
  # 中层循环：遍历每一个组学组合
  for (this_ft in names(ft_combos)) {
    
    ft_list <- ft_combos[[this_ft]] 
    if(this_ft=="Baseline"){
      model_data <- curr_data %>% select(y, all_of(ft_list))
    }else{
      model_data <- curr_data %>% select(y, all_of(cov_names), all_of(ft_list))
    }
    
    # ==========================================================================
    # 使用全局种子生成每个模型特有的种子（保证可重复性）
    # ==========================================================================
    # 为每个任务-特征组合生成唯一的种子
    task_index <- which(names(tasks) == task_name)
    ft_index <- which(names(ft_combos) == this_ft)
    target_seed <- GLOBAL_SEED + task_index * 10000 + ft_index * 100
    
    cat(sprintf('--> Type: %s | Features: %d | Using Generated Seed: %d (from GLOBAL: %d)\n', 
                this_ft, length(ft_list), target_seed, GLOBAL_SEED))
    
    set.seed(target_seed)
    
    fitControl <- trainControl(
      method = "repeatedcv", number = 5, repeats = 5,
      classProbs = TRUE, summaryFunction = twoClassSummary,
      savePredictions = 'all', allowParallel = T, verboseIter = F, sampling = "up"
    )
    
    fit <- train(
      as.formula("y ~ ."), data = model_data, method = "rf",
      trControl = fitControl, tuneLength = 3, metric = 'ROC'
    )
    
    best_mtry <- fit$bestTune$mtry
    
    # 1. 计算 AUC Summary
    temp_roc_calc <- fit$pred %>%
      filter(mtry == best_mtry) %>%
      group_by(Resample) %>%
      summarise(auc = as.numeric(auc(roc(obs, x1, direction='<', quiet=T))), .groups='drop')
    
    current_mean_auc <- mean(temp_roc_calc$auc)
    current_sd_auc   <- sd(temp_roc_calc$auc)
    
    # 保存 Summary
    this_summary <- data.frame(
      outcome = task_name,
      ft_type = this_ft, 
      feature_count = length(ft_list),
      global_seed = GLOBAL_SEED,  # 记录全局种子
      model_seed = target_seed,    # 记录该模型使用的具体种子
      roc_auc = current_mean_auc,
      auc_sd = current_sd_auc
    )
    final_auc_summary <- rbind(final_auc_summary, this_summary)
    
    # 2. 计算 ROC 曲线数据
    best_preds <- fit$pred %>% filter(mtry == best_mtry)
    preds_split <- split(best_preds, best_preds$Resample)
    
    sens_list <- lapply(preds_split, function(df) {
      roc_obj <- pROC::roc(df$obs, df$x1, direction='<', quiet=TRUE)
      vals <- as.numeric(unlist(pROC::coords(
        roc_obj, x = target_spec_grid, input = "specificity", ret = "sensitivity", transpose = FALSE
      )))
      return(vals)
    })
    
    sens_matrix <- do.call(rbind, sens_list)
    mean_sens <- colMeans(sens_matrix, na.rm = TRUE)
    
    mean_roc_interpolated <- data.frame(
      fpr = target_fpr_grid,
      mean_sens = mean_sens
    )
    
    this_roc_data <- cbind(
      outcome = task_name,
      ft_type = this_ft, 
      global_seed = GLOBAL_SEED,
      model_seed = target_seed,
      mean_roc_interpolated
    )
    final_roc_tab <- rbind(final_roc_tab, this_roc_data)
    
    # -------------------------------------------------------------------------
    # 3. 提取变量重要性
    # -------------------------------------------------------------------------
    imp_obj <- varImp(fit, scale = TRUE) 
    
    this_imp_df <- imp_obj$importance
    this_imp_df$Feature <- rownames(this_imp_df)
    rownames(this_imp_df) <- NULL
    
    this_imp_df <- this_imp_df %>%
      mutate(
        outcome = task_name,
        ft_type = this_ft,
        global_seed = GLOBAL_SEED,
        model_seed = target_seed
      ) %>%
      select(outcome, ft_type, Feature, Overall, global_seed, model_seed) %>%
      arrange(desc(Overall))
    
    final_imp_tab <- rbind(final_imp_tab, this_imp_df)
    
    cat(sprintf('   Done. AUC: %.4f (SD: %.4f)\n', current_mean_auc, current_sd_auc))
    
    # 4. 保存模型（如果需要）
    if (task_name == 'PE_vs_NP' && this_ft == 'Multi-Omics') {
      saveRDS(fit, paste0('Anti_results/Best_PE_MultiOmics_Model_GlobalSeed', GLOBAL_SEED, '_ModelSeed', target_seed, '.rds'))
    }
  }
}

stopCluster(cl)

# ==============================================================================
# 保存结果
# ==============================================================================
write.xlsx(final_auc_summary, 'Anti_results/prediction_summary_auc_upsampling_best_0810.xlsx')
write.xlsx(final_roc_tab, 'Anti_results/prediction_roc_curve_data_upsampling_best_0810.xlsx')
write.xlsx(final_imp_tab, 'Anti_results/prediction_importance_upsampling_best.xlsx')




rm(list = ls())
# 请修改为你的实际路径
setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code') 

library(pheatmap); library(openxlsx); library(tidyr); library(RColorBrewer)
library(dplyr); library(ggplot2); library(stringr); library(readxl)
library(Boruta); library(viridis); library(colorspace); library(randomForest)
library(caret); library(pROC); library(reshape2); library(mRMRe)
library(PRROC); library(patchwork); library(doParallel)
source('help_func.R')
source('../R_code/myUtils.R')

# ==============================================================================
# 1. 数据读取与预处理
# ==============================================================================
f_met <- read_excel('Anti_results/Met_sig_merge.xlsx')[['meta_name']]
f_its <- read_excel('Anti_results/ITS_sig_only.xlsx', sheet = 1)[['meta_name']]
f_bac <- read_excel('Anti_results/Taxonomy_s_sig_only.xlsx', sheet = 1)[['meta_name']]

met_dat <- read.csv('data/Weighted_ALL_Met_516.csv', check.names = F) %>% filter(period == 'V1')
its_dat <- read.csv('data/Weighted_ALL_ITS.csv', check.names = F) %>% filter(period == 'V1')
bac_dat <- read.csv('data/Weighted_ALL_Taxonomy.csv', check.names = F) %>% filter(period == 'V1')

met_dat <- filter(met_dat, !is.na(wt2))
met_dat <- met_dat %>% 
  filter(!is.na(HDP)) %>% 
  mutate(
    group = case_when(
      HDP == 1 & preeclampsia == 1 ~ "PE",
      HDP == 1 & preeclampsia == 0 ~ "GH",
      HDP == 0                     ~ "NP",
      TRUE                         ~ NA_character_ 
    )
  )

met_dat$id <- as.numeric(met_dat$id)
its_dat$id <- as.numeric(its_dat$id)
bac_dat$id <- as.numeric(bac_dat$id)

cor_list = c("id", "period", "age", "parity", "edu", "smk", "drk", "wk", "BMI_prep", "group")
cov_data <- met_dat %>% select(all_of(cor_list))
cov_names <- setdiff(names(cov_data), c('id', 'group', 'period'))

data_all <- cov_data %>%
  left_join(met_dat %>% select(id, all_of(f_met)), by = 'id') %>%
  left_join(its_dat %>% select(id, all_of(f_its)), by = 'id') %>%
  left_join(bac_dat %>% select(id, all_of(f_bac)), by = 'id') %>%
  tidyr::drop_na() %>% 
  mutate(
    smk = as.factor(smk),
    drk = as.factor(drk),
    edu = as.factor(edu),
    parity = as.factor(parity),
    group = as.factor(group)
  )

# ==============================================================================
# 2. 外部验证数据准备
# ==============================================================================
valid_raw <- read_excel('../Anti_data/WeBirth_RF_DATA/Multi_omics_only.xlsx')
top10_list <- read_excel('Anti_results/Met_sig_merge.xlsx')[['meta_name']]

cols_in_range <- valid_raw %>% select(MEDP1749:MEDN0502) %>% names()
target_metabolites <- intersect(cols_in_range, top10_list)

valid_raw <- valid_raw %>% select(!all_of(cols_in_range), all_of(target_metabolites))

valid_data_base <- valid_raw %>%
  mutate(smk = 1 - smk) %>%
  mutate(
    smk = as.factor(smk), drk = as.factor(drk), edu = as.factor(edu),
    parity = as.factor(parity), group = as.factor(group)
  )

# ==============================================================================
# 3. 参数配置
# ==============================================================================
# ==============================================================================
# 3. 参数配置 - 修改为全局种子
# ==============================================================================
# 设置全局种子（整个run使用同一个种子）
GLOBAL_SEED <- 123  # 你可以修改这个数字

ft_combos <- list(
  Baseline = cov_names,
  `Met+Fun` = c(f_met, f_its),
  `Met+Bac` = c(f_met, f_bac),
  `Multi-Omics` = c(f_met, f_its, f_bac)
)

tasks <- list(
  PE_vs_NP = list(filter = c('NP', 'PE'), case = 'PE', ctrl = 'NP'),
  Combined_vs_NP = list(filter = c('NP', 'GH', 'PE'), case = c('GH', 'PE'), ctrl = 'NP')
)

# 设置全局随机种子
set.seed(GLOBAL_SEED)

n_outer_repeats <- 5  # 外层重复次数
n_inner_repeats <- 10  # 内层软投票次数

common_1_spec <- seq(0, 1, length.out = 1001)

final_auc_summary <- data.frame() 
final_roc_data_for_plot <- data.frame()

cl <- makeCluster(detectCores() - 2) 
registerDoParallel(cl)

# 采纳轻便、提速的普通 5折 cv
train_control_config <- trainControl(
  method = "cv", number = 5, 
  classProbs = TRUE, summaryFunction = twoClassSummary,
  savePredictions = 'none', allowParallel = TRUE, verboseIter = FALSE
)

# ==============================================================================
# 4. 主循环 - 修改种子生成逻辑
# ==============================================================================
for (task_name in names(tasks)) {
  
  task <- tasks[[task_name]]
  cat(sprintf('\n#### Processing Task: %s ####\n', task_name))
  
  curr_data <- data_all %>% filter(group %in% task$filter)
  curr_data$y <- factor(ifelse(curr_data$group %in% task$case, 'x1', 'x0'), levels = c('x1', 'x0'))
  
  curr_valid <- valid_data_base %>% filter(group %in% task$filter)
  curr_valid$y <- factor(ifelse(curr_valid$group %in% task$case, 'x1', 'x0'), levels = c('x1', 'x0'))
  curr_valid <- curr_valid %>% tidyr::drop_na()
  
  for (ft_type in names(ft_combos)) {
    
    ft_list <- ft_combos[[ft_type]]
    if(ft_type=="Baseline"){
      features_needed <- c(ft_list)
    } else {
      features_needed <- c(cov_names, ft_list)
    }
    
    model_data <- curr_data %>% select(y, all_of(features_needed))
    
    valid_ready <- TRUE
    if (!all(features_needed %in% colnames(curr_valid))) {
      valid_ready <- FALSE
      cat(sprintf("  [Warning] Missing features in validation set for %s. Skipping...\n", ft_type))
      next
    }
    
    # --- 参数网格定义 ---
    relaxed_group <- c("Fungi", "Bacteria", "Fun+Bac")
    if (ft_type == "Baseline") {
      current_grid_params <- expand.grid(maxnodes = 100, nodesize = 2, ntree = 500, mtry_val = 3)
    } else if(ft_type %in% relaxed_group){
      current_grid_params <- expand.grid(maxnodes = 50, nodesize = 3, ntree = 500, mtry_val = 5)
    } else {
      current_grid_params <- expand.grid(maxnodes = 25, nodesize = 10, ntree = 500, mtry_val = 6)
    }
    
    mn <- current_grid_params$maxnodes[1]
    ns <- current_grid_params$nodesize[1]
    nt <- current_grid_params$ntree[1]
    my_tune_grid <- data.frame(mtry = current_grid_params$mtry_val[1])
    
    cat(sprintf('\n--> Type: %s | Feats: %d | Using GLOBAL seed: %d\n', 
                ft_type, length(features_needed), GLOBAL_SEED))
    
    outer_auc_values <- c()
    outer_sens_matrix <- matrix(NA, nrow = n_outer_repeats, ncol = length(common_1_spec))
    
    # ==========================================================================
    # Outer Loop: 外层实验 - 使用全局种子生成可重复的序列
    # ==========================================================================
    for (run_idx in 1:n_outer_repeats) {
      
      # 使用全局种子生成不同的子种子（保证可重复性）
      main_seed <- GLOBAL_SEED + run_idx * 1000 + which(names(tasks) == task_name) * 10000 + which(names(ft_combos) == ft_type) * 100000
      accum_valid_probs <- rep(0, nrow(curr_valid))
      
      # --- Inner Loop: Soft Voting (内层子模型) ---
      for (iter in 1:n_inner_repeats) {
        current_sub_seed <- main_seed + iter
        set.seed(current_sub_seed)
        
        fit_model <- train(
          as.formula("y ~ ."), data = model_data, method = "rf",
          trControl = train_control_config, tuneGrid = my_tune_grid, metric = 'ROC',
          ntree = nt, maxnodes = mn, nodesize = ns
        )
        
        preds_ext <- predict(fit_model, newdata = curr_valid, type = "prob")[, "x1"]
        accum_valid_probs <- accum_valid_probs + preds_ext
      }
      
      # 计算单次 Outer Run 的集成预测结果
      mean_valid_probs <- accum_valid_probs / n_inner_repeats
      roc_ext <- pROC::roc(curr_valid$y, mean_valid_probs, levels = c('x0', 'x1'), direction = '<', quiet = TRUE)
      
      # 收集 AUC
      current_run_auc <- as.numeric(pROC::auc(roc_ext))
      outer_auc_values <- c(outer_auc_values, current_run_auc)
      
      # 收集曲线插值结果，存入矩阵用于最后求均值
      coords_res <- pROC::coords(roc = roc_ext, x = "all", ret = c("specificity", "sensitivity"), transpose = FALSE)
      interp_res <- approx(x = 1 - coords_res$specificity, y = coords_res$sensitivity, xout = common_1_spec, ties = mean)
      outer_sens_matrix[run_idx, ] <- interp_res$y
      
      cat(sprintf('    [Run %d] Ext AUC: %.4f\n', run_idx, current_run_auc))
    } # 结束 Outer Run
    
    # -------------------------------------------------------------------------
    # 汇总计算 Mean AUC / SD
    # -------------------------------------------------------------------------
    final_mean_ext_auc <- mean(outer_auc_values)
    final_sd_ext_auc   <- sd(outer_auc_values)
    
    cat(sprintf('  >> Finished! Mean Ext AUC: %.4f (SD: %.4f)\n', final_mean_ext_auc, final_sd_ext_auc))
    
    # 计算外层实验的 Mean 曲线
    final_mean_sens <- colMeans(outer_sens_matrix, na.rm = TRUE)
    
    # 补充 (0,0) 坐标
    tmp_plot_df <- data.frame(
      outcome = task_name,
      ft_type = ft_type,
      global_seed = GLOBAL_SEED,
      fpr = c(0, common_1_spec),
      sens = c(0, final_mean_sens)
    ) %>% distinct()
    
    final_roc_data_for_plot <- rbind(final_roc_data_for_plot, tmp_plot_df)
    
    summary_row <- data.frame(
      outcome = task_name,
      ft_type = ft_type,
      n_experiments = n_outer_repeats,
      global_seed = GLOBAL_SEED,
      mean_auc = final_mean_ext_auc,
      sd_auc = final_sd_ext_auc
    )
    final_auc_summary <- rbind(final_auc_summary, summary_row)
    
  } # End Feature Type Loop
} # End Task Loop

stopCluster(cl)

# ==============================================================================
# 5. 结果保存
# ==============================================================================
write.xlsx(final_auc_summary, 'Anti_results/Final_Summary_Stats_WeBirth_0810.xlsx')
write.csv(final_roc_data_for_plot, 'Anti_results/Final_ROC_Curves_Data_WeBirth_0810.csv', row.names = FALSE)

