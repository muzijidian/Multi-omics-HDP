library(dplyr)
library(ggplot2)
library(openxlsx)
library(networkdata)
library(igraph)
library(ggraph)
library(graphlayouts)
library(ggforce)


focus_16S_info <- read.xlsx('./Fig3b/16S_sig.xlsx')  %>%  rename(micro_name=name)
focus_16S_info$meta_name <- paste0('bac_', focus_16S_info$meta_name)
focus_16S_info$source <- '16S'

focus_ITS_info <- read.xlsx('./data/ITS_sig.xlsx')  %>%  rename(micro_name=name)
focus_ITS_info$meta_name <- paste0('fun_', focus_ITS_info$meta_name)
focus_ITS_info$source <- 'ITS'

focus_MGS_info <- read.xlsx('./data/Taxonomy_s_sig.xlsx')  %>%  rename(micro_name=name)
focus_MGS_info$meta_name <- paste0('tax_', focus_MGS_info$meta_name)
focus_MGS_info$source <- 'Metagenomic species'

THSBC_16S_ITS <- read.xlsx('./data/HDP_Subgroup_Pearson_ALL_ITSV1_16SV1_c_wt1_antibiotic_ascov.xlsx')
THSBC_MGS_ITS <- read.xlsx('./data/HDP_Subgroup_Pearson_ALL_ITSV1_TaxsV1_c_wt2_antibiotic_ascov.xlsx')
THSBC_16S_MGS <- read.xlsx('./data/HDP_Subgroup_Pearson_ALL_16SV1_TaxsV1_c_wt2_antibiotic_ascov.xlsx')

THSBC_16S_ITS <- THSBC_16S_ITS  %>%  left_join(focus_16S_info, by=c('Var1'='meta_name'))  %>%  left_join(focus_ITS_info,  by=c('Var2'='meta_name')) 
THSBC_MGS_ITS <- THSBC_MGS_ITS  %>%  left_join(focus_ITS_info, by=c('Var1'='meta_name'))  %>%  left_join(focus_MGS_info,  by=c('Var2'='meta_name'))  
THSBC_16S_MGS <- THSBC_16S_MGS  %>%  left_join(focus_16S_info, by=c('Var1'='meta_name'))  %>%  left_join(focus_MGS_info,  by=c('Var2'='meta_name'))

THSBC_16S_ITS <- na.omit(THSBC_16S_ITS)
THSBC_MGS_ITS <- na.omit(THSBC_MGS_ITS)
THSBC_16S_MGS <- na.omit(THSBC_16S_MGS)
THSBC_Met_Micros <- rbind(THSBC_16S_ITS, THSBC_MGS_ITS, THSBC_16S_MGS)
edges_all <- THSBC_Met_Micros  %>%  select(x=micro_name.x, y=micro_name.y, beta=r, p=pv)
edges_all$p_FDR <- p.adjust(edges_all$p, method = 'fdr')

name <- c(THSBC_Met_Micros$micro_name.x, THSBC_Met_Micros$micro_name.y)
group <- c(THSBC_Met_Micros$source.x, THSBC_Met_Micros$source.y)
nodes_df <- data.frame(name=name, group=group)
nodes_df <- unique(nodes_df)
nodes_df <- nodes_df  %>%  group_by(group) %>% arrange(group, name) %>% ungroup()

# Build igraph object
edges_all_filtered <- edges_all  %>%  filter(p<0.05)
network <- graph_from_data_frame(
    d = edges_all_filtered,
    vertices = nodes_df,
    directed = FALSE
)

# Calculate nodes degree
node_degrees <- degree(network)
V(network)$degree <- node_degrees

# Set base position for node groups
base_positions <- list(
    "Metagenomic species" = c(0, 6),    
    "16S" = c(-3, -6),    
    "ITS" = c(4, -2)      
)

n_points <- table(nodes_df$group)
make_layout_matrix <- function(meta){
    angles <- seq(0, 2 * pi, length.out = n_points[[meta]] + 1)[-1]
    layout_coords <- data.frame(
        x = base_positions[[meta]][1] + n_points[[meta]] / 4 * cos(angles),
        y = base_positions[[meta]][2] + n_points[[meta]] / 4 * sin(angles)
    )
    return(layout_coords)
}

layout_matrix <- rbind(make_layout_matrix('16S'), make_layout_matrix('ITS'),make_layout_matrix('Metagenomic species'))

p <- ggraph(network, layout = "manual", x = layout_matrix[,1], y = layout_matrix[,2]) +
    geom_edge_link0(aes(
        alpha = ifelse(p<0.05, 0.5, 0.3),
        width = abs(beta),
        color = ifelse(beta<0, 'blue', 'red')
    )) +
    geom_node_point(aes(
        fill = group,
        size = degree
    ), shape = 21, color = "white", stroke = 0.5) +
    geom_node_text(
        aes(label = name), 
        repel = TRUE, 
        size = 2,    
        fontface = "bold",  
        bg.color = "white",  
        bg.r = 0.1,         
        segment.color = "grey50",  
        max.overlaps = Inf,  
        show.legend = FALSE) +
    geom_mark_hull(
        aes(x = x, y = y, fill = group, group = group),
        concavity = 2,
        expand = unit(4, "mm"),
        alpha = 0.15
    ) +
    scale_edge_width_continuous(range = c(0.5, 2.5)) +  
    scale_edge_alpha_continuous(range = c(0.3, 0.5)) +  
    scale_edge_color_manual(
        values = c("blue" = "#4575B4", "red" = "#D73027", 'gray'='gray'),
        labels = c("Negative", "Positive"),
        name = "Correlation"
    ) +
    scale_fill_manual(
        values = c("16S" = "#7FC97F", 
                  "ITS" = "#BEAED4", 
                  "Metagenomic species" = "#FB9A99"
                  )  
    ) +
    theme_graph(base_family = "serif") +
    coord_fixed(clip = "off") +
    xlim(c(-12, 12)) +
    ylim(c(-12, 12)) +
    labs(
        edge_alpha = "Significance (-log10 p-value)",
        edge_width = "Correlation strength",
        fill = "Group"
    )

# print(p)
