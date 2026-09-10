rm(list = ls())
library(dplyr)
library(ggplot2)
library(openxlsx)
library(stringr)
library(colorspace)
library(ggbreak)  
library(ggnewscale) # [New] Used to map a second set of colors in the same plot (bottom category color blocks)
setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code')

# 1. Read and merge data, and perform FDR correction ---------------------------------------
df_bac = read.xlsx('Anti_results/variance_explain_bacteria_adonis2_longitudinal_new.xlsx') %>% 
  mutate(Microbiome = "Bacteria",
         padj = p.adjust(pv_margin, method = "fdr"))

df_fun = read.xlsx('Anti_results/variance_explain_fungi_adonis2_longitudinal_new.xlsx') %>% 
  mutate(Microbiome = "Fungi",
         padj = p.adjust(pv_margin, method = "fdr"))

df = bind_rows(df_bac, df_fun)

# 2. Variable renaming and significance labeling ------------------------------------------------
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
                                'drk' ~ "Alcohol consumption",
                                'tpa_score5' ~ "Physical activity", 
                                'sleep_score' ~ "Sleep score",
                                
                                # --- Medication use ---
                                'aspirin_painkiller_use' ~ "Aspirin/analgesic use",
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
         
         # [Modification 1] Use FDR-corrected padj to generate significance asterisks
         text = case_when(padj <= 0.001 ~ '***',
                          padj <= 0.01 ~ '**',
                          padj <= 0.05 ~ '*', 
                          .default = '')
  )

# 3. Calculate global ranking (sorted in descending order by the sum of R2 of Bacteria and Fungi) ------------------

# [Modification 2] Discard the original fixed order, automatically calculate the X-axis factor order in descending order of R2
var_order <- df_clean %>%
  group_by(var_label) %>%
  summarise(total_R2 = sum(R2_margin, na.rm = TRUE)) %>%
  arrange(desc(total_R2)) %>%
  pull(var_label)

df_clean$var_label <- factor(df_clean$var_label, levels = var_order)

# Lock Category and Microbiome factors
df_clean$Category = factor(df_clean$Category, 
                           levels = c("Maternal baseline characteristics", 
                                      "Lifestyle factors", 
                                      "Medication use", 
                                      "Clinical history", 
                                      "Dietary intake"))
df_clean$Microbiome <- factor(df_clean$Microbiome, levels = c("Bacteria", "Fungi"))

# 4. Plotting --------------------------------------------------------------

# Custom colors (bar chart)
color_bac <- colorspace::lighten('#FFB2B9', amount = 0.1) 
color_fun <- colorspace::lighten('#ebdcb2', amount = 0.1) 

# Custom colors (bottom Category color blocks, using low-saturation macaron colors to avoid overshadowing)
cat_colors <- c("Maternal baseline characteristics" = "#e9ac70",
                "Lifestyle factors" = "#e0756e",
                "Medication use" = "#687f99",
                "Clinical history" = "#9fc07f",
                "Dietary intake" = "#9c91b8")



p1 <- ggplot(df_clean, aes(x = var_label)) +
  
  # [Modification] Widen the bars: increase width from 0.8 to 0.9, and adjust position_dodge accordingly
  geom_col(aes(y = R2_margin * 100, fill = Microbiome), 
           position = position_dodge(width = 0.9), width = 0.9, alpha = 0.9) +
  
  # [Modification] Adjust asterisk position to match the widened bars
  geom_text(aes(y = R2_margin * 100 - 0.01, label = text, group = Microbiome), 
            position = position_dodge(width = 0.9), vjust = 0, size = 4.2, color = "black", 
            fontface = "bold") +
  
  # [Modification] Manually set new colors, and set the legend to a single column and order it first
  scale_fill_manual(values = c("Bacteria" = color_bac, "Fungi" = color_fun),
                    guide = guide_legend(ncol = 1, order = 1)) +
  
  # ==========================================
# Introduce a second set of fill mappings
ggnewscale::new_scale_fill() +
  
  # [Modification] Narrow the Category blocks: decrease height from 0.05 to 0.03
  geom_tile(aes(y = -0.05, fill = Category), height = 0.03, width = 0.9) +
  
  # [Modification] Set Category legend to a single column and order it second
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
    
    # [Modification] Ensure both legends are displayed horizontally side-by-side at the top
    # legend.position = c(0.5, 0.8),
    # legend.position = "right",
    legend.box = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = "none",
    plot.margin = unit(c(5, 5, 5, 20), "mm")
  )

# The subsequent y-axis break and save code remains unchanged
p2 <- p1 + scale_y_break(breaks = c(0.37, 3.5), scales = 0.5)
ggsave('Anti_graph/fig2f_variance_explained_Combined_Sorted_Break_R2_margin.pdf', p2, width = 9, height = 4.7)








rm(list = ls())

library(dplyr)
library(ggplot2)
library(openxlsx)
library(stringr)
library(colorspace)
# [Modification] Removed ggbreak package
library(ggnewscale) 

# =========================================================================
# 1. Global parameters and color settings
# =========================================================================

# Custom colors for T1, T2, T3 (using soft macaron color palette)
period_colors <- c("T1" = "#D9C59A",
                   "T2" = "#A8D8B9",  # Soft green
                   "T3" = "#7EABCB"   # Soft blue
)  # Soft yellow (T1)

# Custom colors for bottom Category color blocks

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
# 2. Define the integrated [Data processing + Plotting] function
# =========================================================================

plot_variance_by_period <- function(file_path, title_name, output_pdf) {
  
  # --- 1. Read data and perform FDR correction ---
  df <- read.xlsx(file_path) %>% 
    # Perform FDR correction grouped by period
    group_by(Period) %>%
    mutate(padj = p.adjust(pv_margin, method = "BH")) %>%
    ungroup()
  
  # --- 2. Variable renaming and significance labeling ---
  df_clean <- df %>% 
    mutate(var_label = case_match(var_name, 
                                  'wk' ~ "Gestational age",
                                  'age' ~ "Maternal age",
                                  'BMI_prep' ~ "Pre-pregnancy BMI",
                                  'parity' ~ "Parity",
                                  'edu' ~ "Education level", 
                                  'smk' ~ "Smoking status",
                                  'drk' ~ "Alcohol consumption",
                                  'tpa_score5' ~ "Physical activity", 
                                  'sleep_score' ~ "Sleep score",
                                  'aspirin_painkiller_use' ~ "Aspirin/analgesic use",
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
  
  # --- 3. Calculate ranking and fix factors ---
  # Sort in descending order based on the sum of R2_margin across T1, T2, T3
  var_order <- df_clean %>%
    group_by(var_label) %>%
    summarise(total_R2 = sum(R2_margin, na.rm = TRUE)) %>%
    arrange(desc(total_R2)) %>%
    pull(var_label)
  
  df_clean$var_label <- factor(df_clean$var_label, levels = var_order)
  df_clean$Category <- factor(df_clean$Category, levels = cat_levels)
  df_clean$Period <- factor(df_clean$Period, levels = c("T1", "T2", "T3"))
  
  # --- 4. Plotting ---
  p1 <- ggplot(df_clean, aes(x = var_label)) +
    
    # Bar chart: grouped side-by-side by period
    geom_col(aes(y = R2_margin * 100, fill = Period), 
             position = position_dodge(width = 0.9), width = 0.9, alpha = 0.9) +
    
    # Significance asterisks: aligned with the bars of each period
    geom_text(aes(y = R2_margin * 100 - 0.01, label = text, group = Period), 
              position = position_dodge(width = 0.9), vjust = 0, size = 7, color = "black") +
    
    # First set of legends: T1, T2, T3
    scale_fill_manual(values = period_colors,
                      guide = guide_legend(ncol = 1, order = 1, title = "Trimester")) +
    
    # === Introduce the second set of legends (bottom Category) ===
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
  
  # --- 5. Save and return the plot ---
  # [Modification] Do not generate p2 anymore, directly save and return p1
  ggsave(output_pdf, p1, width = 14, height = 4.5)
  
  return(p1)
}

# =========================================================================
# 3. Execute plotting and output PDF
# =========================================================================

# Plot for Fungi
p_fungi <- plot_variance_by_period(
  file_path = 'Anti_results/variance_explain_fungi_adonis2_T123.xlsx', 
  title_name = "Fungi", 
  output_pdf = 'Anti_graph/fig2f_variance_explained_Fungi_T123.pdf'
)

# Plot for Bacteria
p_bacteria <- plot_variance_by_period(
  file_path = 'Anti_results/variance_explain_bacteria_adonis2_T123.xlsx', 
  title_name = "Bacteria", 
  output_pdf = 'Anti_graph/fig2f_variance_explained_Bacteria_T123.pdf'
)

print(p_fungi)
print(p_bacteria)