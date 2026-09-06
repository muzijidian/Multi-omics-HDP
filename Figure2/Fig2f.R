rm(list = ls())

# 1. 加载依赖包 -----------------------------------------------------------
library(dplyr)
library(data.table)
library(stringr)
library(reshape2)
library(openxlsx)
library(ggplot2)
library(vegan)
library(lme4)
library(lmerTest)
library(doParallel)
library(patchwork)
library(haven) 

# 设置工作目录
setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code')
select = dplyr::select

# 2. 数据读取与预处理 --------------------------------------------------------

# --- A. 读取菌群数据 ---
# bcdist = read.csv('data/Weighted_ALL_Taxonomy.csv') # 根据需要切换
bcdist = read.csv('data/Weighted_ALL_ITS.csv')

bcdist <- bcdist[!is.na(bcdist$wt1), ]
bcdist$period <- sub("^V", "T", as.character(bcdist$period))

bcdist <- bcdist %>% 
  filter(!is.na(HDP)) %>% 
  mutate(
    group = case_when(
      HDP == 1 & preeclampsia == 1 ~ "PE",
      HDP == 1 & preeclampsia == 0 ~ "GH",
      HDP == 0                    ~ "NP",
      TRUE                        ~ NA_character_
    )
  )

load('data/beta_diversity_rawdata_forplot.RData')
# bacteia_list = bact_mgs_s_uni
bacteia_list = fung_its_uni

# --- B. 读取并合并附加数据 ---

# 1. Disease Data (静态)
diease_data = read.xlsx('../data/THSBC_new_data.xlsx')
bcdist$id <- as.character(bcdist$id)
diease_data$id <- as.character(diease_data$id)

disease_subset <- diease_data %>%
  select(id, aspirin_painkiller_use, no_antibiotic_use, 
         hpt_prep, diabetes_prep, gdm_ever, ghpt_ever, preeclam_ever, family_diabetes, family_hpt)

# 2. Lifestyle Data (静态/新引入)
data_lifestyle <- read_dta("data/lifestyle factor.dta")
data_lifestyle$id <- as.character(data_lifestyle$id)
lifestyle_subset <- data_lifestyle %>% 
  select(id, tpa_score5, sleep_score)

# 3. Dietary Data (纵向 - 移除 _a，使用全孕期数据)
data_dietary <- readRDS("Anti_results/longtab_diet_cov_7281.rds")
data_dietary$id <- as.character(data_dietary$id)

if("period" %in% names(data_dietary)){
  data_dietary$period <- sub("^V", "T", as.character(data_dietary$period))
} else {
  stop("Error: data_dietary 中缺少 period 列")
}

data_dietary_subset <- data_dietary %>%
  select(id, period, cat_egg, cat_carbo, cat_meat, cat_dairy, cat_vege, cat_fruit)

# 4. 合并数据
cat("合并前 bcdist 行数:", nrow(bcdist), "\n")

datac <- bcdist %>% 
  left_join(data_dietary_subset, by = c("id", "period")) %>% # 纵向双键匹配 (此时的饮食数据也是动态的)
  left_join(disease_subset, by = "id") %>%                   # 静态单键匹配
  left_join(lifestyle_subset, by = "id")                     # 新增：静态单键匹配

cat("合并后 datac 行数:", nrow(datac), "\n")

# --- C. 定义所有变量及其分类 (重新构建类别) ---

# 提取出 smk 和 drk，原类别更名为 Maternal baseline characteristics
baseline_list   = c('age', 'parity', 'edu', 'wk', 'BMI_prep') 
# 新增 Lifestyle 类别
lifestyle_list  = c('smk', 'drk', 'tpa_score5', 'sleep_score')
medication_list = c('aspirin_painkiller_use', 'no_antibiotic_use')
clinical_list   = c('diabetes_prep', 'gdm_ever', 'ghpt_ever', 'preeclam_ever', "family_diabetes", 'family_hpt')
# 饮食变量去掉 _a
dietary_list    = c('cat_egg', 'cat_carbo', 'cat_meat', 'cat_dairy', 'cat_vege', 'cat_fruit')

# 构建变量映射表
var_map <- rbind(
  data.frame(var = baseline_list,   category = "Maternal baseline characteristics"),
  data.frame(var = lifestyle_list,  category = "Lifestyle factors"),
  data.frame(var = medication_list, category = "Medication use"),
  data.frame(var = clinical_list,   category = "Clinical history"),
  data.frame(var = dietary_list,    category = "Dietary intake")
)

# --- D. 最终筛选与清洗 ---

keep_cols <- c("id", "period", "group", 
               var_map$var, 
               all_of(bacteia_list))

datac <- datac[, keep_cols, drop = FALSE]
datac = na.omit(datac) 
cat("最终纳入分析的样本量 (NA removed): ", nrow(datac), "\n")


# 4. 核心分析循环 (只看 All 人群，不再分组) -------------------------------------------

set.seed(1111)
result_combined = data.frame() 

for (t in c('T1', 'T2', 'T3')){
  
  # 数据切片
  if (t != 'All'){
    data_t = datac %>% filter(period==t)
  } else{
    data_t = datac
  }
  
  cat('\n================ Processing Period:', t, 'Samples:', nrow(data_t), '================\n')
  
  # 准备距离矩阵 (只计算全人群，去掉了 PE/GH/NP)
  vdist_data = data_t %>% select(all_of(bacteia_list))
  vdist = vegdist(vdist_data, 'bray')
  
  # 循环所有单变量
  for (i in 1:nrow(var_map)){
    
    cov_name = var_map$var[i]
    cov_cat  = var_map$category[i]
    
    model_terms = cov_name
    if (cov_name=='wk' & t=='All') {
      model_terms = c('wk','period')
    }
    
    formula_str = paste0('vdist ~ ', paste(model_terms, collapse = "+"))
    
    # 运行 adonis2
    # 针对 All 模式加入了 strata=data_t$id 以控制重复测量带来的影响
    if (t == 'All') {
      fit_all = adonis2(as.formula(formula_str), data = data_t, 
                        permutations = 999, parallel = 8)
    } else {
      fit_all = adonis2(as.formula(formula_str), data = data_t, 
                        permutations = 999, parallel = 8)
    }
    
    # 存储结果 (精简版：仅保留 All)
    result_combined = rbind(result_combined, data.frame(
      period = t,
      Category = cov_cat, 
      var_name = cov_name,
      type = "Single",
      R2_All = fit_all$R2[1], 
      pv_All = fit_all$`Pr(>F)`[1]
    ))
  }
}

# 5. 结果展示与保存 ---------------------------------------------------------

cat("\nAnalysis Completed.\n")
print(head(result_combined))

write.xlsx(result_combined, 'Anti_results/variance_explain_fungi_T123.xlsx')








rm(list = ls())
library(dplyr)
library(ggplot2)
library(openxlsx)
library(stringr)
library(colorspace)
library(ggbreak)  
library(ggnewscale) # 【新增】用于在同一张图里映射第二套颜色（底部分类色块）
setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code')

# 1. 读取并合并数据，并执行 FDR 矫正 ---------------------------------------
df_bac = read.xlsx('Anti_results/variance_explain_bacteria_adonis2_longitudinal_new.xlsx') %>% 
  mutate(Microbiome = "Bacteria",
         padj = p.adjust(pv_margin, method = "fdr"))

df_fun = read.xlsx('Anti_results/variance_explain_fungi_adonis2_longitudinal_new.xlsx') %>% 
  mutate(Microbiome = "Fungi",
         padj = p.adjust(pv_margin, method = "fdr"))

df = bind_rows(df_bac, df_fun)

# 2. 变量重命名与显著性标记 ------------------------------------------------
df_clean = df %>% 
  mutate(var_label = case_match(var_name, 
                                # --- Maternal baseline characteristics ---
                                'wk' ~ "Gestational age",
                                'age' ~ "Maternal age",
                                'BMI_prep' ~ "Pre-pregnancy BMI",
                                'parity' ~ "Parity",
                                'edu' ~ "Education level", 
                                
                                # --- Lifestyle factors ---
                                'smk' ~ "Smoking status",
                                'drk' ~ "Drinking status",
                                'tpa_score5' ~ "Physical activity", 
                                'sleep_score' ~ "Sleep score",
                                
                                # --- Medication use ---
                                'aspirin_painkiller_use' ~ "Aspirin/Painkiller",
                                'no_antibiotic_use' ~ "Antibiotic use",
                                
                                # --- Clinical history ---
                                'gdm_ever' ~ "History of GDM", 
                                'ghpt_ever' ~ "History of HDP", 
                                'family_diabetes' ~ "Family Diabetes", 
                                'family_hpt' ~ "Family Hypertension",
                                
                                # --- Dietary intake ---
                                'cat_vege' ~ "Vegetables", 
                                'cat_fruit' ~ "Fruits",
                                'cat_meat' ~ "Meat",
                                'cat_egg' ~ "Eggs", 
                                'cat_dairy' ~ "Dairy", 
                                'cat_carbo' ~ "Total fast carbohydrates", 
                                
                                .default = var_name),
         
         # 【修改1】改用经过 FDR 矫正的 padj 来生成显著性星号
         text = case_when(padj <= 0.001 ~ '***',
                          padj <= 0.01 ~ '**',
                          padj <= 0.05 ~ '*', 
                          .default = '')
  )

# 3. 计算全局排序 (按 Bacteria 和 Fungi 的 R2 总和降序排列) ------------------

# 【修改2】抛弃原有的定死顺序，自动按 R2 降序计算 X 轴因子顺序
var_order <- df_clean %>%
  group_by(var_label) %>%
  summarise(total_R2 = sum(R2_margin, na.rm = TRUE)) %>%
  arrange(desc(total_R2)) %>%
  pull(var_label)

df_clean$var_label <- factor(df_clean$var_label, levels = var_order)

# 锁定 Category 和 Microbiome 因子
df_clean$Category = factor(df_clean$Category, 
                           levels = c("Maternal baseline characteristics", 
                                      "Lifestyle factors", 
                                      "Medication use", 
                                      "Clinical history", 
                                      "Dietary intake"))
df_clean$Microbiome <- factor(df_clean$Microbiome, levels = c("Bacteria", "Fungi"))

# 4. 绘图 --------------------------------------------------------------

# 自定义颜色 (柱状图)
color_bac <- colorspace::lighten('#FFB2B9', amount = 0.1) 
color_fun <- colorspace::lighten('#ebdcb2', amount = 0.1) 

# 自定义颜色 (底部 Category 色块，采用低饱和度马卡龙色以防喧宾夺主)
cat_colors <- c("Maternal baseline characteristics" = "#e9ac70",
                "Lifestyle factors" = "#e0756e",
                "Medication use" = "#687f99",
                "Clinical history" = "#9fc07f",
                "Dietary intake" = "#9c91b8")



p1 <- ggplot(df_clean, aes(x = var_label)) +
  
  # 【修改】条形柱变宽：width 从 0.8 增加到 0.9，并相应调整 position_dodge
  geom_col(aes(y = R2_margin * 100, fill = Microbiome), 
           position = position_dodge(width = 0.9), width = 0.9, alpha = 0.9) +
  
  # 【修改】调整星号位置以匹配加宽后的柱子
  geom_text(aes(y = R2_margin * 100 - 0.01, label = text, group = Microbiome), 
            position = position_dodge(width = 0.9), vjust = 0, size = 4.2, color = "black", 
            fontface = "bold") +
  
  # 【修改】手动设置新颜色，并设置图例为单列且排在第一
  scale_fill_manual(values = c("Bacteria" = color_bac, "Fungi" = color_fun),
                    guide = guide_legend(ncol = 1, order = 1)) +
  
  # ==========================================
# 引入第二套 fill 映射
ggnewscale::new_scale_fill() +
  
  # 【修改】Category 的方块变窄：height 从 0.05 减小到 0.03
  geom_tile(aes(y = -0.05, fill = Category), height = 0.03, width = 0.9) +
  
  # 【修改】设置 Category 图例为单列且排在第二
  scale_fill_manual(values = cat_colors,
                    guide = guide_legend(ncol = 1, order = 2)) +
  # ==========================================

labs(x = '', y = 'Variance explained (%)') +
  theme_classic(base_line_size = 0.3) +
  theme(
    axis.text.y = element_text(color = 'black', size = 14),
    axis.title = element_text(size = 15),
    axis.line = element_line(linewidth = 0.5, color = "black"), 
    axis.text.x = element_text(size = 13, color = 'black', angle = 90, hjust = 1, 
                               margin = ggplot2::margin(t = 5)),
    
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.line.x.top = element_blank(),
    
    # 【修改】确保两个图例水平并排显示在顶部
    # legend.position = c(0.5, 0.8),
    # legend.position = "right",
    legend.box = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = "none",
    plot.margin = unit(c(5, 5, 5, 20), "mm")
  )

# 后续的 y 轴截断和保存代码保持不变
p2 <- p1 + scale_y_break(breaks = c(0.37, 3.5), scales = 0.5)
ggsave('Anti_graph/fig2f_variance_explained_Combined_Sorted_Break_R2_margin.pdf', p2, width = 9, height = 4.7)








rm(list = ls())

library(dplyr)
library(ggplot2)
library(openxlsx)
library(stringr)
library(colorspace)
# 【修改】已移除 ggbreak 包
library(ggnewscale) 

# =========================================================================
# 1. 全局参数与颜色设置
# =========================================================================

# 自定义 T1, T2, T3 的颜色 (使用柔和的马卡龙色系)
period_colors <- c("T1" = "#D9C59A",
                   "T2" = "#A8D8B9",  # 柔和绿
                   "T3" = "#7EABCB"  # 柔和蓝
)  # 柔和黄

# 自定义底部 Category 色块颜色

cat_colors <- c("Maternal baseline characteristics" = "#e9ac70",
                "Lifestyle factors" = "#e0756e",
                "Medication use" = "#687f99",
                "Clinical history" = "#9fc07f",
                "Dietary intake" = "#9c91b8")

cat_levels <- c("Maternal baseline characteristics", 
                "Lifestyle factors", 
                "Medication use", 
                "Clinical history", 
                "Dietary intake")

# =========================================================================
# 2. 定义【数据处理+绘图】一体化函数
# =========================================================================

plot_variance_by_period <- function(file_path, title_name, output_pdf) {
  
  # --- 1. 读取数据与 FDR 矫正 ---
  df <- read.xlsx(file_path) %>% 
    # 按照 period 分组进行 FDR 矫正
    group_by(Period) %>%
    mutate(padj = p.adjust(pv_margin, method = "BH")) %>%
    ungroup()
  
  # --- 2. 变量重命名与显著性标记 ---
  df_clean <- df %>% 
    mutate(var_label = case_match(var_name, 
                                  'wk' ~ "Gestational age",
                                  'age' ~ "Maternal age",
                                  'BMI_prep' ~ "Pre-pregnancy BMI",
                                  'parity' ~ "Parity",
                                  'edu' ~ "Education level", 
                                  'smk' ~ "Smoking status",
                                  'drk' ~ "Drinking status",
                                  'tpa_score5' ~ "Physical activity", 
                                  'sleep_score' ~ "Sleep score",
                                  'aspirin_painkiller_use' ~ "Aspirin/Painkiller",
                                  'no_antibiotic_use' ~ "Antibiotic use",
                                  'gdm_ever' ~ "History of GDM", 
                                  'ghpt_ever' ~ "History of HDP", 
                                  'family_diabetes' ~ "Family diabetes", 
                                  'family_hpt' ~ "Family hypertension",
                                  'cat_vege' ~ "Vegetables", 
                                  'cat_fruit' ~ "Fruits",
                                  'cat_meat' ~ "Meat",
                                  'cat_egg' ~ "Eggs", 
                                  'cat_dairy' ~ "Dairy", 
                                  'cat_carbo' ~ "Total fast carbohydrates", 
                                  .default = var_name),
           
           text = case_when(padj <= 0.001 ~ '***',
                            padj <= 0.01 ~ '**',
                            padj <= 0.05 ~ '*', 
                            .default = ''))
  
  # --- 3. 计算排序并固定因子 ---
  # 按该变量在 T1, T2, T3 的 R2_All 总和进行降序排列
  var_order <- df_clean %>%
    group_by(var_label) %>%
    summarise(total_R2 = sum(R2_margin, na.rm = TRUE)) %>%
    arrange(desc(total_R2)) %>%
    pull(var_label)
  
  df_clean$var_label <- factor(df_clean$var_label, levels = var_order)
  df_clean$Category <- factor(df_clean$Category, levels = cat_levels)
  df_clean$Period <- factor(df_clean$Period, levels = c("T1", "T2", "T3"))
  
  # --- 4. 绘图 ---
  p1 <- ggplot(df_clean, aes(x = var_label)) +
    
    # 柱状图：按 period 分组并列
    geom_col(aes(y = R2_margin * 100, fill = Period), 
             position = position_dodge(width = 0.9), width = 0.9, alpha = 0.9) +
    
    # 显著性星号：按 period 的柱子对齐
    geom_text(aes(y = R2_margin * 100 - 0.01, label = text, group = Period), 
              position = position_dodge(width = 0.9), vjust = 0, size = 7, color = "black") +
    
    # 第一套图例：T1, T2, T3
    scale_fill_manual(values = period_colors,
                      guide = guide_legend(ncol = 1, order = 1, title = "Trimester")) +
    
    # === 引入第二套图例 (底部 Category) ===
    ggnewscale::new_scale_fill() +
    geom_tile(aes(y = -0.05, fill = Category), height = 0.03, width = 0.9) +
    scale_fill_manual(values = cat_colors,
                      guide = guide_legend(ncol = 1, order = 2)) +
    # ======================================
  
  labs(x = '', y = 'Variance explained (%)', title = title_name) +
    theme_classic(base_line_size = 0.3) +
    theme(
      # plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.title = element_blank(),
      axis.text.y = element_text(color = 'black', size = 15),
      axis.title = element_text(size = 16),
      axis.line = element_line(linewidth = 0.5, color = "black"), 
      axis.text.x = element_text(size = 14, color = 'black', angle = 35, hjust = 1, 
                                 margin = ggplot2::margin(t = 5)),
      
      axis.text.x.top = element_blank(),
      axis.ticks.x.top = element_blank(),
      axis.line.x.top = element_blank(),
      
      # legend.position = "top", 
      # legend.box = "horizontal",
      # legend.text = element_text(size = 12),
      # legend.title = element_text(size = 13, face = "bold"),
      legend.position = "none", 
      plot.margin = unit(c(5, 5, 5, 20), "mm")
    )
  
  # --- 5. 保存并返回图形 ---
  # 【修改】不再生成 p2，直接保存并返回 p1
  ggsave(output_pdf, p1, width = 14, height = 4.5)
  
  return(p1)
}

# =========================================================================
# 3. 执行绘图并输出 PDF
# =========================================================================

# 画真菌 (Fungi) 的图
p_fungi <- plot_variance_by_period(
  file_path = 'Anti_results/variance_explain_fungi_adonis2_T123.xlsx', 
  title_name = "Fungi", 
  output_pdf = 'Anti_graph/fig2f_variance_explained_Fungi_T123.pdf'
)

# 画细菌 (Bacteria) 的图
p_bacteria <- plot_variance_by_period(
  file_path = 'Anti_results/variance_explain_bacteria_adonis2_T123.xlsx', 
  title_name = "Bacteria", 
  output_pdf = 'Anti_graph/fig2f_variance_explained_Bacteria_T123.pdf'
)

print(p_fungi)
print(p_bacteria)



