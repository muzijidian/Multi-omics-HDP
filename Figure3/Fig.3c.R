library(dplyr)
library(stringr)
library(reshape2)
library(openxlsx)
library(ggplot2)
z_trans <- function(x){
  (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
}

# load data -----------------------------------
ft_16S = read.xlsx('../data/16S_sig_info.xlsx')
ft_16S = ft_16S %>% select(data_source, source, ft_name = meta_name, beta, pv_adj_FDR, text = name)
ft_ITS = read.xlsx('../data/ITS_sig_info.xlsx')
ft_ITS = ft_ITS %>% select(data_source, source, ft_name = meta_name, beta, pv_adj_FDR, text = name)
ft_mgen_g = read.xlsx('../data/Taxonomy_g_sig_info.xlsx')
ft_mgen_g = ft_mgen_g %>% select(source, ft_name = meta_name, beta, pv_adj_FDR, text = name) %>% mutate(data_source = 'mg_g')
ft_mgen_s = read.xlsx('../data/Taxonomy_s_sig_info.xlsx')
ft_mgen_s = ft_mgen_s %>% select(source, ft_name = meta_name, beta, pv_adj_FDR, text = name) %>% mutate(data_source = 'mg_s')

ft_diff = bind_rows(ft_16S, ft_ITS, ft_mgen_g, ft_mgen_s)
ft_diff = ft_diff %>% mutate(ft_name = paste0(ft_name, '_', data_source)) %>% 
  mutate(outcome = str_remove(source, '_T\\d+'),
         period = str_extract(source, 'T\\d+'))

T1_list = union(ft_diff %>% filter(period == 'T1' & pv_adj_FDR < 0.1) %>% pull(ft_name),
                ft_diff %>% filter(data_source=='ITS' & period == 'T2' & pv_adj_FDR < 0.1) %>% pull(ft_name))
ft_diff = ft_diff %>% mutate(type = ifelse(ft_name %in% T1_list, 'T1', 'T3'))

# load abundance data -----------------------------
value_16S = read.xlsx('../data/HDP_16S_arcsinZ_prevalence_0.25_minabund_1e-4_BPinfo.xlsx')
value_16S = value_16S %>% select(id:GH, matches('^g\\d+')) %>% 
  reshape2::melt(id.vars = 1:6, variable.name = 'ft_name') %>% 
  mutate(ft_name = paste0(ft_name, '_16S'))
mean_16S = value_16S %>% group_by(ft_name, period) %>% 
  summarise(value = weighted.mean(value, 1/prob)) %>% 
  left_join(value_16S %>% filter(HDP==1) %>% 
              group_by(ft_name) %>% 
              mutate(value = z_trans(value)) %>% 
              group_by(ft_name, period) %>% 
              summarise(value_HDP = weighted.mean(value, 1/prob))) %>% 
  left_join(value_16S %>% filter(HDP==0) %>% 
              group_by(ft_name) %>% 
              mutate(value = z_trans(value)) %>% 
              group_by(ft_name, period) %>% 
              summarise(value_nHDP = weighted.mean(value, 1/prob)))

value_ITS = read.xlsx('../data/HDP_ITS_arcsinZ_prevalence_0.25_minabund_1e-4_BPinfo.xlsx')
value_ITS = value_ITS %>% select(id:GH, matches('^g\\d+')) %>% 
  reshape2::melt(id.vars = 1:6, variable.name = 'ft_name') %>% 
  mutate(ft_name = paste0(ft_name, '_ITS'))
mean_ITS = value_ITS %>% group_by(ft_name, period) %>% 
  summarise(value = weighted.mean(value, 1/prob)) %>% 
  left_join(value_ITS %>% filter(HDP==1) %>% 
              group_by(ft_name) %>% 
              mutate(value = z_trans(value)) %>% 
              group_by(ft_name, period) %>% 
              summarise(value_HDP = weighted.mean(value, 1/prob))) %>% 
  left_join(value_ITS %>% filter(HDP==0) %>% 
              group_by(ft_name) %>% 
              mutate(value = z_trans(value)) %>% 
              group_by(ft_name, period) %>% 
              summarise(value_nHDP = weighted.mean(value, 1/prob)))

value_mgen_g = read.xlsx('../data/HDP_Taxonomy_g_arcsinZ_prevalence_0.25_minabund_1e-4_BPinfo.xlsx')
value_mgen_g = value_mgen_g %>% select(id:GH, matches('^g\\d+')) %>% 
  reshape2::melt(id.vars = 1:6, variable.name = 'ft_name') %>% 
  mutate(ft_name = paste0(ft_name, '_mg_g'))
mean_mgen_g = value_mgen_g %>% group_by(ft_name, period) %>% 
  summarise(value = weighted.mean(value, 1/prob)) %>% 
  left_join(value_mgen_g %>% filter(HDP==1) %>% 
              group_by(ft_name) %>% 
              mutate(value = z_trans(value)) %>% 
              group_by(ft_name, period) %>% 
              summarise(value_HDP = weighted.mean(value, 1/prob))) %>% 
  left_join(value_mgen_g %>% filter(HDP==0) %>% 
              group_by(ft_name) %>% 
              mutate(value = z_trans(value)) %>% 
              group_by(ft_name, period) %>% 
              summarise(value_nHDP = weighted.mean(value, 1/prob)))

value_mgen_s = read.xlsx('../data/HDP_Taxonomy_s_arcsinZ_prevalence_0.25_minabund_1e-4_BPinfo.xlsx')
value_mgen_s = value_mgen_s %>% select(id:GH, matches('^s\\d+')) %>% 
  reshape2::melt(id.vars = 1:6, variable.name = 'ft_name') %>% 
  mutate(ft_name = paste0(ft_name, '_mg_s'))
mean_mgen_s = value_mgen_s %>% group_by(ft_name, period) %>% 
  summarise(value = weighted.mean(value, 1/prob)) %>% 
  left_join(value_mgen_s %>% filter(HDP==1) %>% 
              group_by(ft_name) %>% 
              mutate(value = z_trans(value)) %>% 
              group_by(ft_name, period) %>% 
              summarise(value_HDP = weighted.mean(value, 1/prob))) %>% 
  left_join(value_mgen_s %>% filter(HDP==0) %>% 
              group_by(ft_name) %>% 
              mutate(value = z_trans(value)) %>% 
              group_by(ft_name, period) %>% 
              summarise(value_nHDP = weighted.mean(value, 1/prob)))

ft_mean_all = bind_rows(mean_16S, mean_ITS, mean_mgen_g, mean_mgen_s) %>% 
  ungroup() %>% 
  mutate(diff = ifelse(ft_name %in% unique(ft_diff$ft_name), 1, 0)) %>% 
  mutate(diff2 = case_when(ft_name %in% T1_list ~ '1', 
                           diff == 1 ~'2', .default = '0'))
# table(ft_mean_all$diff2)/3

# trend cluster all ------------------------------------
library(latrend)

ft_mean_all$ID = as.character(as.numeric(factor(ft_mean_all$ft_name)))
trend_data = ft_mean_all %>% 
  mutate(time = as.numeric(str_remove(period, 'V')))
  # select(-value) %>% rename(value = value)
# length(unique(trend_data$ID)) # 383

options(latrend.id = "ID", latrend.time = "time")
# plotTrajectories(trend_data, response = "value")
set.seed(1234)
kmlMethod = lcMethodKML(response = "value", nClusters = 3, nbRedrawing = 20)
kmlMethods = lcMethods(kmlMethod, nClusters = 1:12)
model = latrendBatch(kmlMethods, data = trend_data)
# 最佳聚类数
plotMetric(model, c("ASW", "BIC", "WMAE", "WRSS"))+
  geom_vline(xintercept = 4, linetype=2, color='coral')+
  theme_bw()
ggsave('graph/kml_cluster_bestK.pdf', width = 4.4, height = 4)

best_fit = subset(model, nClusters == 4, drop = TRUE)
plot(best_fit)+
  scale_x_continuous(expand = c(0,0.2), breaks = 1:3, labels = c('T1','T2','T3'))
cluster_traj = clusterTrajectories(best_fit)
cluster = data.frame(ID = ids(best_fit), Cluster = trajectoryAssignments(best_fit))
df = trend_data %>% left_join(cluster)

df_unique = df %>% select(ft_name, diff, Cluster) %>% unique()
chisq.test(df_unique$diff, df_unique$Cluster)
count = table(df_unique$diff, df_unique$Cluster)
ratio = count[2,]/ count[1,]
df_text = data.frame(Cluster = colnames(count), 
                     count0 = count[1,], count1 = count[2,],
                     ratio = ratio)

cluster_label = paste0(paste0('Cluster ',names(table(df$Cluster))), 
                       ' (n=',table(df$Cluster)/3,')')
names(cluster_label) = names(table(df$Cluster))
# 差异变化 cluster聚类图
class_col = c('#1e996a','#3F5D7D','#E6B91E','#e96e40')
names(class_col) = names(cluster_label)
ggplot(df, aes(x = time, y = value))+
  facet_wrap(.~Cluster, labeller = labeller(Cluster=cluster_label))+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.5, color='grey80')+
  geom_line(aes(group = ID, color = Cluster), alpha = 0.04)+
  geom_line(aes(group = Cluster, color = Cluster), data=cluster_traj, linewidth=0.8)+
  geom_text(aes(label = sprintf('Differential: %d, Neutral: %d', count1, count0)), 
            x = 3.1, y = 0.61, data = df_text, hjust = 1, size = 4)+
  geom_text(aes(label = sprintf('Ratio: %.2f', ratio)), 
            x = 3.1, y = 0.4, data = df_text, hjust = 1, size = 4)+
  scale_color_manual(values = class_col)+
  scale_x_continuous(expand = c(0,0.2), breaks = 1:3, labels = c('T1','T2','T3'))+
  labs(y = 'Relative abundance (arcsine Z)')+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10.2, color = "black", 
                                    vjust = 1.9, hjust = 0.5, angle = 90),
        axis.text.x = element_text(size = 10.2, color = "black"), 
        axis.text.y = element_text(size = 9.2, color = "grey10", 
                                   vjust = 0.5, hjust = 0.5, angle = 0),
        axis.ticks.length = unit(1.5, 'mm'),
        strip.text = element_text(size = 11, color = 'black'),
        strip.background = element_rect(fill = 'grey90'),
        legend.position = 'none')
ggsave('graph/micro_abundance_cluster.pdf', width = 5, height = 5*2/3, device = cairo_pdf)

################ not run ####################
# trend cluster T1-------------------------------
trend_data = ft_diff %>% filter(type == 'T1') %>% 
  rename(ID = ft_name, value = beta) %>% 
  mutate(time = as.numeric(factor(source)))
# length(unique(trend_data$ID)) # 18

options(latrend.id = "ID", latrend.time = "time")
plotTrajectories(trend_data, response = "value")
# kml 选择cluster数
set.seed(1234)
kmlMethod = lcMethodKML(response = "value", nClusters = 3, nbRedrawing = 20)
kmlMethods = lcMethods(kmlMethod, nClusters = 1:17)
model = latrendBatch(kmlMethods, data = trend_data)
# 最佳聚类数
plotMetric(model, c("Dunn", "ASW", "BIC", "WMAE","WRSS", "estimationTime"))+
  # geom_vline(xintercept = 11, linetype=2, color='coral')+
  # scale_x_continuous(breaks = seq(0,25,5))+
  theme_bw()
# ggsave('graph/kml_cluster_bestK.png', width = 8, height = 5)

best_fit = subset(model, nClusters == 17, drop = TRUE)
plot(best_fit) # 重新排序
cluster_traj = clusterTrajectories(best_fit) %>% 
  left_join(unique(trend_data %>% select(time, source))) %>% 
  mutate(outcome = str_remove(source, '_T\\d+'),
         period = str_extract(source, 'T\\d+'))
cluster = data.frame(ID = ids(best_fit), Cluster = trajectoryAssignments(best_fit))
df = trend_data %>% left_join(cluster)

cluster_label = paste0(paste0('Cluster ',names(table(df$Cluster))), 
                       ' (n=',table(df$Cluster)/9,')')
names(cluster_label) = names(table(df$Cluster))
# 差异变化 cluster聚类图
set.seed(9999)
class_col = sample(cols4all::c4a('rainbow_wh_rd', length(cluster_label)))
names(class_col) = names(cluster_label)
ggplot(df, aes(x = period, y = value))+
  facet_grid(outcome~Cluster, labeller = labeller(Cluster=cluster_label))+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.5, color='grey80')+
  geom_line(aes(group = ID, color = Cluster), alpha = 0.15)+
  geom_line(aes(group = Cluster, color = Cluster), data=cluster_traj, linewidth=0.45)+
  scale_color_manual(values = class_col)+
  scale_x_discrete(expand = c(0,0.2))+
  labs(y = 'log OR')+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10.2, color = "black", 
                                    vjust = 1.9, hjust = 0.5, angle = 90),
        axis.text.x = element_text(size = 10.2, color = "black"), 
        axis.text.y = element_text(size = 9.2, color = "grey10", 
                                   vjust = 0.5, hjust = 0.5, angle = 0),
        axis.ticks.length = unit(1.5, 'mm'),
        strip.text = element_text(size = 11, color = 'black'),
        strip.background = element_rect(fill = 'grey90'),
        legend.position = 'none')
# ggsave('graph/fig3_diff_trend.pdf', width = 4.5, height = 3)

# trend cluster T3-------------------------------
trend_data = ft_diff %>% filter(type == 'T3') %>% 
  rename(ID = ft_name, value = beta) %>% 
  mutate(time = as.numeric(factor(source)))
# length(unique(trend_data$ID)) # 64

options(latrend.id = "ID", latrend.time = "time")
# plotTrajectories(trend_data, response = "value")
# kml 选择cluster数
set.seed(1234)
kmlMethod = lcMethodKML(response = "value", nClusters = 3, nbRedrawing = 20)
kmlMethods = lcMethods(kmlMethod, nClusters = 1:25)
model = latrendBatch(kmlMethods, data = trend_data)
# 最佳聚类数
plotMetric(model, c("Dunn", "ASW", "BIC", "WMAE","WRSS", "estimationTime"))+
  geom_vline(xintercept = 7, linetype=2, color='coral')+
  scale_x_continuous(breaks = seq(0,25,5))+
  theme_bw()
# ggsave('graph/kml_cluster_bestK_T3.png', width = 8, height = 5)

best_fit = subset(model, nClusters == 7, drop = TRUE)
# plot(best_fit) # 重新排序
cluster_traj = clusterTrajectories(best_fit) %>% 
  left_join(unique(trend_data %>% select(time, source))) %>% 
  mutate(outcome = str_remove(source, '_T\\d+'),
         period = str_extract(source, 'T\\d+'))
cluster = data.frame(ID = ids(best_fit), Cluster = trajectoryAssignments(best_fit))
df = trend_data %>% left_join(cluster) %>% 
  mutate(outcome = str_remove(source, '_T\\d+'),
         period = str_extract(source, 'T\\d+'))
vi = df %>% group_by(Cluster, time) %>% 
  summarise(mean_value = mean(value)) %>% 
  left_join(cluster_traj)


cluster_label = paste0(paste0('Cluster ',names(table(df$Cluster))), 
                       ' (n=',table(df$Cluster)/9,')')
names(cluster_label) = names(table(df$Cluster))


set.seed(1000)
class_col = sample(cols4all::c4a('rainbow_wh_rd', length(cluster_label)))
names(class_col) = names(cluster_label)
ggplot(df, aes(x = period, y = value))+
  facet_grid(outcome~Cluster, labeller = labeller(Cluster=cluster_label))+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.5, color='grey80')+
  geom_line(aes(group = ID, color = Cluster), alpha = 0.1)+
  geom_line(aes(group = Cluster, color = Cluster), data=cluster_traj, linewidth=0.5)+
  scale_color_manual(values = class_col)+
  scale_x_discrete(expand = c(0,0.2))+
  labs(y = 'log OR')+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10.2, color = "black", 
                                    vjust = 1.9, hjust = 0.5, angle = 90),
        axis.text.x = element_text(size = 10.2, color = "black"), 
        axis.text.y = element_text(size = 9.2, color = "grey10", 
                                   vjust = 0.5, hjust = 0.5, angle = 0),
        axis.ticks.length = unit(1.5, 'mm'),
        strip.text = element_text(size = 11, color = 'black'),
        strip.background = element_rect(fill = 'grey90'),
        legend.position = 'none')
ggsave('graph/fig3_diff_trend.pdf', width = 4.5, height = 3)