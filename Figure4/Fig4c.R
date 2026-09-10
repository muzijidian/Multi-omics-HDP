rm(list = ls())
library(ggplot2)
library(openxlsx)
library(dplyr)
library(ggpubr)
setwd("/Users/mzjd/Documents/HDP-multiomics/haonan_code")
set.seed(123)

colorset <- c("#FFFFB3", "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", 
              "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#3498DB", "#FEB29B")
class_order <- c("AlAm", "AAD", "BA", "BSD", "CHD", "FA", "GP", "HRC", "NUC", "OAD", "SPL", "Others")

# Create class color mapping
names(colorset) <- class_order

hdp_t1 <- read.xlsx('../Anti_data/preprocessed/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov_class3.xlsx', sheet = "HDP_T1")
pe_t1 <- read.xlsx('../Anti_data/preprocessed/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov_class3.xlsx', sheet = "PE_T1")
gh_t1 <- read.xlsx('../Anti_data/preprocessed/HDP_ALL_Met_wt1_batch_Info_no_antibiotic_as_cov_class3.xlsx', sheet = "GH_T2")

# Filter sugg_cmpd_name with pv_adj_FDR < 0.05
hdp_t1_sel <- hdp_t1 %>% filter(pv_adj_FDR < 0.05) %>% select(sugg_cmpd_name, beta, Class.I)
pe_t1_sel <- pe_t1 %>% filter(pv_adj_FDR < 0.05) %>% select(sugg_cmpd_name, beta, Class.I)
gh_t1_sel <- gh_t1 %>% filter(pv_adj_FDR < 0.05) %>% select(sugg_cmpd_name, beta, Class.I)

# Read HBC (GAOLING) and WeBirth (Westlake) data
gaoling_t1 <- read.xlsx('../Anti_data/preprocessed/HBC_HDP_ALL_Met_Info_no_antibiotic_as_cov.xlsx', sheet = "HDP_T1")
westlake_t1 <- read.xlsx('../Anti_data/preprocessed/WeBirth_HDP_ALL_Met_Info_no_antibiotic_as_cov.xlsx', sheet = "HDP_T1")

gaoling_pe <- read.xlsx('../Anti_data/preprocessed/HBC_HDP_ALL_Met_Info_no_antibiotic_as_cov.xlsx', sheet = "PE_T1")
westlake_pe <- read.xlsx('../Anti_data/preprocessed/WeBirth_HDP_ALL_Met_Info_no_antibiotic_as_cov.xlsx', sheet = "PE_T1")

gaoling_gh <- read.xlsx('../Anti_data/preprocessed/HBC_HDP_ALL_Met_Info_no_antibiotic_as_cov.xlsx', sheet = "GH_T1")
westlake_gh <- read.xlsx('../Anti_data/preprocessed/WeBirth_HDP_ALL_Met_Info_no_antibiotic_as_cov.xlsx', sheet = "GH_T1")

# Intersection of HDP and HBC (GAOLING)
hdp_t1_gaoling_common <- hdp_t1_sel %>%
  filter(sugg_cmpd_name %in% gaoling_t1$sugg_cmpd_name) %>%
  left_join(gaoling_t1 %>% select(sugg_cmpd_name, beta_gaoling = beta), by = "sugg_cmpd_name")

# Intersection of HDP and WeBirth (WESTLAKE)  
hdp_t1_westlake_common <- hdp_t1_sel %>%
  filter(sugg_cmpd_name %in% westlake_t1$sugg_cmpd_name) %>%
  left_join(westlake_t1 %>% select(sugg_cmpd_name, beta_westlake = beta), by = "sugg_cmpd_name")

# Intersection of PE and HBC (GAOLING)
pe_t1_gaoling_common <- pe_t1_sel %>%
  filter(sugg_cmpd_name %in% gaoling_pe$sugg_cmpd_name) %>%
  left_join(gaoling_pe %>% select(sugg_cmpd_name, beta_gaoling = beta), by = "sugg_cmpd_name")

# Intersection of PE and WeBirth (WESTLAKE)
pe_t1_westlake_common <- pe_t1_sel %>%
  filter(sugg_cmpd_name %in% westlake_pe$sugg_cmpd_name) %>%
  left_join(westlake_pe %>% select(sugg_cmpd_name, beta_westlake = beta), by = "sugg_cmpd_name")

# Intersection of GH and HBC (GAOLING)
gh_t1_gaoling_common <- gh_t1_sel %>%
  filter(sugg_cmpd_name %in% gaoling_gh$sugg_cmpd_name) %>%
  left_join(gaoling_gh %>% select(sugg_cmpd_name, beta_gaoling = beta), by = "sugg_cmpd_name")

# Intersection of GH and WeBirth (WESTLAKE)
gh_t1_westlake_common <- gh_t1_sel %>%
  filter(sugg_cmpd_name %in% westlake_gh$sugg_cmpd_name) %>%
  left_join(westlake_gh %>% select(sugg_cmpd_name, beta_westlake = beta), by = "sugg_cmpd_name")

hdp_t1_gaoling <- hdp_t1_gaoling_common %>% 
  select(sugg_cmpd_name, beta, beta_source = beta_gaoling, Class.I) %>% 
  mutate(Cohorts = "HBC", period = "HDP_T1")

hdp_t1_westlake <- hdp_t1_westlake_common %>% 
  select(sugg_cmpd_name, beta, beta_source = beta_westlake, Class.I) %>% 
  mutate(Cohorts = "WeBirth", period = "HDP_T1")

pe_t1_gaoling <- pe_t1_gaoling_common %>% 
  select(sugg_cmpd_name, beta, beta_source = beta_gaoling, Class.I) %>% 
  mutate(Cohorts = "HBC", period = "PE_T1")

pe_t1_westlake <- pe_t1_westlake_common %>% 
  select(sugg_cmpd_name, beta, beta_source = beta_westlake, Class.I) %>% 
  mutate(Cohorts = "WeBirth", period = "PE_T1")

gh_t1_gaoling <- gh_t1_gaoling_common %>% 
  select(sugg_cmpd_name, beta, beta_source = beta_gaoling, Class.I) %>% 
  mutate(Cohorts = "HBC", period = "GH_T1")

gh_t1_westlake <- gh_t1_westlake_common %>% 
  select(sugg_cmpd_name, beta, beta_source = beta_westlake, Class.I) %>% 
  mutate(Cohorts = "WeBirth", period = "GH_T1")

beta_comb <- rbind(
  pe_t1_gaoling,
  pe_t1_westlake,
  hdp_t1_gaoling,
  hdp_t1_westlake
  # gh_t1_gaoling,
  # gh_t1_westlake
)

# Modify the order of the period column
beta_comb$period <- factor(beta_comb$period, levels = c('PE_T1', 'HDP_T1', 'GH_T1'))
beta_comb$Class.I <- factor(beta_comb$Class.I, levels = class_order)

actual_classes <- unique(beta_comb$Class.I)
colorset_filtered <- colorset[names(colorset) %in% actual_classes]

cohort_shapes <- c("HBC" = 16, "WeBirth" = 17)
period_col <- c("HBC" = "#2E5266", "WeBirth" = "#F44336")

p1 <- ggplot(beta_comb, aes(x = beta, y = beta_source, color = Class.I, shape = Cohorts)) +
  geom_point(alpha = 1, size = 2.5) +
  # Add regression line for HBC cohort
  geom_smooth(data = subset(beta_comb, Cohorts == "HBC"), 
              method = 'lm', alpha = 0.15, color = period_col["HBC"], 
              fill = period_col["HBC"], se = TRUE) +
  # Add regression line for WeBirth cohort
  geom_smooth(data = subset(beta_comb, Cohorts == "WeBirth"), 
              method = 'lm', alpha = 0.15, color = period_col["WeBirth"], 
              fill = period_col["WeBirth"], se = TRUE) +
  scale_color_manual(values = colorset_filtered, name = "Metabolite Class") +
  scale_shape_manual(values = cohort_shapes, name = "Cohorts") +
  facet_wrap(~period, nrow = 1, scales = 'free') +
  stat_cor(data = subset(beta_comb, Cohorts == "HBC"),
           method = "pearson", r.accuracy = 0.01, p.accuracy = 0.001,
           label.x.npc = 0.2, label.y.npc = 0.13, 
           show.legend = FALSE, color = period_col["HBC"], size = 5) +
  # Add correlation statistics for WeBirth cohort, using WeBirth colors
  stat_cor(data = subset(beta_comb, Cohorts == "WeBirth"),
           method = "pearson", r.accuracy = 0.01, p.accuracy = 0.001,
           label.x.npc = 0.2, label.y.npc = 0.23, 
           show.legend = FALSE, color = period_col["WeBirth"], size = 5) +
  labs(x = 'Beta coefficients in THSBC', 
       y = 'Beta coefficients in WeBirth/HBC') +
  theme_classic() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(0.99), hjust = 0, face = 'bold'),
        axis.title = element_text(size = rel(1.3)),
        axis.text = element_text(size = rel(1.3),  color = 'black'),
        legend.position = "bottom",
        legend.box = "horizontal") +
  guides(
    color = guide_legend(
      override.aes = list(size = 3), 
      ncol = length(actual_classes),  # Keep Metabolite Classes arranged horizontally
      title = "Metabolite Classes"
    ),
    shape = "none"
  ) +
  theme(
    legend.box = "vertical",
    legend.box.just = "left",
    legend.box.spacing = unit(0.15, "cm"),     # Reduce spacing between the two legend groups
    legend.spacing.y = unit(0.05, "cm"),       # Reduce vertical spacing between legend items
    legend.spacing.x = unit(0.1, "cm"),        # Reduce horizontal spacing between legend items
    legend.margin = ggplot2::margin(2, 2, 2, 2, "pt"),  # Reduce outer margin
    legend.text = element_text(size = 9),       # Shrink text size
    legend.title = element_text(size = 10, face = "plain")  # Shrink title size
  )

# Save plot
ggsave('Anti_graph/Fig4c_pearson_crossdata_class_HDP_PE_class3.pdf', width = 4.5, height = 4)