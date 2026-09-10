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
# Set global seed 
# ==============================================================================
GLOBAL_SEED <- 555 
set.seed(GLOBAL_SEED)

# --- Load feature lists ---
f_met <- read_excel('Anti_results/Met_sig.xlsx')[['meta_name']]
f_its <- read_excel('Anti_results/ITS_sig.xlsx', sheet = 1)[['meta_name']]
f_bac <- read_excel('Anti_results/Taxonomy_s_sig.xlsx', sheet = 1)[['meta_name']]

# --- Load data ---
met_dat <- read.csv('data/Weighted_ALL_Met_516.csv', check.names = F) %>% 
  filter(period == 'V1')
its_dat <- read.csv('data/Weighted_ALL_ITS.csv', check.names = F) %>% 
  filter(period == 'V1')
bac_dat <- read.csv('data/Weighted_ALL_Taxonomy.csv', check.names = F) %>% 
  filter(period == 'V1')

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

# --- Initialize result containers ---
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
  
  # Middle loop: Iterate through each omics combination
  for (this_ft in names(ft_combos)) {
    
    ft_list <- ft_combos[[this_ft]] 
    if(this_ft == "Baseline"){
      model_data <- curr_data %>% select(y, all_of(ft_list))
    } else {
      model_data <- curr_data %>% select(y, all_of(cov_names), all_of(ft_list))
    }
    
    task_index <- which(names(tasks) == task_name)
    ft_index <- which(names(ft_combos) == this_ft)
    target_seed <- GLOBAL_SEED + task_index * 10000 + ft_index * 100
    
    cat(sprintf('--> Type: %s | Features: %d\n', this_ft, length(ft_list)))
    
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
    
    # 1. Calculate AUC Summary
    temp_roc_calc <- fit$pred %>%
      filter(mtry == best_mtry) %>%
      group_by(Resample) %>%
      summarise(auc = as.numeric(auc(roc(obs, x1, direction='<', quiet=T))), .groups='drop')
    
    current_mean_auc <- mean(temp_roc_calc$auc)
    current_sd_auc   <- sd(temp_roc_calc$auc)
    
    # Save Summary
    this_summary <- data.frame(
      outcome = task_name,
      ft_type = this_ft, 
      feature_count = length(ft_list),
      roc_auc = current_mean_auc,
      auc_sd = current_sd_auc
    )
    final_auc_summary <- rbind(final_auc_summary, this_summary)
    
    # 2. Calculate ROC curve data
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
      mean_roc_interpolated
    )
    final_roc_tab <- rbind(final_roc_tab, this_roc_data)
    
    # -------------------------------------------------------------------------
    # 3. Extract variable importance
    # -------------------------------------------------------------------------
    imp_obj <- varImp(fit, scale = TRUE) 
    
    this_imp_df <- imp_obj$importance
    this_imp_df$Feature <- rownames(this_imp_df)
    rownames(this_imp_df) <- NULL
    
    this_imp_df <- this_imp_df %>%
      mutate(
        outcome = task_name,
        ft_type = this_ft
      ) %>%
      select(outcome, ft_type, Feature, Overall) %>%
      arrange(desc(Overall))
    
    final_imp_tab <- rbind(final_imp_tab, this_imp_df)
    
    cat(sprintf('   Done. AUC: %.4f (SD: %.4f)\n', current_mean_auc, current_sd_auc))
    
  }
}

stopCluster(cl)

# ==============================================================================
# Save results
# ==============================================================================

write.xlsx(final_auc_summary, 'Anti_results/THSBC_AUC.xlsx')
write.xlsx(final_roc_tab, 'Anti_results/THSBC_ROC.xlsx')
write.xlsx(final_imp_tab, 'Anti_results/THSBC_importance.xlsx')



# ==============================================================================
# SECOND PART OF THE SCRIPT
# ==============================================================================
rm(list = ls())

# Please modify to your actual path
setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code') 

library(pheatmap); library(openxlsx); library(tidyr); library(RColorBrewer)
library(dplyr); library(ggplot2); library(stringr); library(readxl)
library(Boruta); library(viridis); library(colorspace); library(randomForest)
library(caret); library(pROC); library(reshape2); library(mRMRe)
library(PRROC); library(patchwork); library(doParallel)
source('help_func.R')
source('../R_code/myUtils.R')

# ==============================================================================
# 1. Data reading and preprocessing
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
# 2. External validation data preparation
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
# 3. Parameter configuration
# ==============================================================================

# Set global seed
GLOBAL_SEED <- 123 

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

set.seed(GLOBAL_SEED)

n_outer_repeats <- 5  # Outer loop repeats
n_inner_repeats <- 10 # Inner soft voting repeats

common_1_spec <- seq(0, 1, length.out = 1001)

final_auc_summary <- data.frame() 
final_roc_data_for_plot <- data.frame()

cl <- makeCluster(detectCores() - 2) 
registerDoParallel(cl)

# Adopt lightweight, accelerated standard 5-fold CV
train_control_config <- trainControl(
  method = "cv", number = 5, 
  classProbs = TRUE, summaryFunction = twoClassSummary,
  savePredictions = 'none', allowParallel = TRUE, verboseIter = FALSE
)

# ==============================================================================
# 4. Main loop
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
    if(ft_type == "Baseline"){
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
    
    # --- Parameter grid definition ---
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
    
    cat(sprintf('\n--> Type: %s | Feats: %d\n', ft_type, length(features_needed)))
    
    outer_auc_values <- c()
    outer_sens_matrix <- matrix(NA, nrow = n_outer_repeats, ncol = length(common_1_spec))
    
    # ==========================================================================
    # Outer Loop: Outer experiment - Generate reproducible sequences using a global seed
    # ==========================================================================
    for (run_idx in 1:n_outer_repeats) {
      
      # Generate different sub-seeds using the global seed (ensuring reproducibility)
      main_seed <- GLOBAL_SEED + run_idx * 1000 + which(names(tasks) == task_name) * 10000 + which(names(ft_combos) == ft_type) * 100000
      accum_valid_probs <- rep(0, nrow(curr_valid))
      
      # --- Inner Loop: Soft Voting (Inner sub-models) ---
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
      
      # Calculate the ensemble prediction result for a single Outer Run
      mean_valid_probs <- accum_valid_probs / n_inner_repeats
      roc_ext <- pROC::roc(curr_valid$y, mean_valid_probs, levels = c('x0', 'x1'), direction = '<', quiet = TRUE)
      
      # Collect AUC
      current_run_auc <- as.numeric(pROC::auc(roc_ext))
      outer_auc_values <- c(outer_auc_values, current_run_auc)
      
      # Collect curve interpolation results, save to matrix for final averaging
      coords_res <- pROC::coords(roc = roc_ext, x = "all", ret = c("specificity", "sensitivity"), transpose = FALSE)
      interp_res <- approx(x = 1 - coords_res$specificity, y = coords_res$sensitivity, xout = common_1_spec, ties = mean)
      outer_sens_matrix[run_idx, ] <- interp_res$y
      
      cat(sprintf('    [Run %d] Ext AUC: %.4f\n', run_idx, current_run_auc))
    } # End Outer Run
    
    # -------------------------------------------------------------------------
    # Aggregate and calculate Mean AUC / SD
    # -------------------------------------------------------------------------
    final_mean_ext_auc <- mean(outer_auc_values)
    final_sd_ext_auc   <- sd(outer_auc_values)
    
    cat(sprintf('  >> Finished! Mean Ext AUC: %.4f (SD: %.4f)\n', final_mean_ext_auc, final_sd_ext_auc))
    
    # Calculate Mean curve for the outer experiment
    final_mean_sens <- colMeans(outer_sens_matrix, na.rm = TRUE)
    
    # Supplement (0,0) coordinates
    tmp_plot_df <- data.frame(
      outcome = task_name,
      ft_type = ft_type,
      fpr = c(0, common_1_spec),
      sens = c(0, final_mean_sens)
    ) %>% distinct()
    
    final_roc_data_for_plot <- rbind(final_roc_data_for_plot, tmp_plot_df)
    
    summary_row <- data.frame(
      outcome = task_name,
      ft_type = ft_type,
      n_experiments = n_outer_repeats,
      mean_auc = final_mean_ext_auc,
      sd_auc = final_sd_ext_auc
    )
    final_auc_summary <- rbind(final_auc_summary, summary_row)
    
  } # End Feature Type Loop
} # End Task Loop

stopCluster(cl)

# ==============================================================================
# 5. Save results
# ==============================================================================
write.xlsx(final_auc_summary, 'Anti_results/WeBirth_validation_AUC_PE.xlsx')
write.csv(final_roc_data_for_plot, 'Anti_results/WeBirth_validation_ROC_PE.csv', row.names = FALSE)