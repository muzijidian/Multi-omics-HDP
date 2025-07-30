
rm(list = ls())
setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code')

library(grid)
library(forestploter)
library(readxl)

dt = read_excel("../data/G130_forest.xlsx")

# 计算置信区间
dt$CI_lower <- dt$beta - 1.96 * dt$se
dt$CI_upper <- dt$beta + 1.96 * dt$se
dt$` ` <- paste(rep(" ", 15), collapse = " ")


dt$Interaction <- ifelse(!is.na(dt$interaction), 
                         sprintf("%.3f", dt$interaction), 
                         "")


gh_row <- which(dt$Age == "GH (48/18)")
pe_row <- which(dt$Age == "PE (27/11)")


blank_row <- dt[1,]
blank_row[1,] <- NA
blank_row$Age <- ""
blank_row$` ` <- paste(rep(" ", 17), collapse = " ")
blank_row$Interaction <- ""


gh_last_row <- gh_row + 2


dt_new <- rbind(
  dt[1:gh_last_row, ],
  blank_row,
  dt[(gh_last_row+1):nrow(dt), ]
)


pe_row <- pe_row + 1


tm <- forest_theme(base_size = 16,
                   core = list(padding = unit(c(7, 11), "mm")), 
                   refline_gp = gpar(lwd = 2,            # 参考线的线宽
                                     lty = "dashed",     # 参考线的线型
                                     col = "#656c8a"      # 参考线的颜色
                   ),
                   vertline_lwd = 2,         # 垂直线的线宽，可以添加一条额外的垂直线，如果没有就不显示
                   vertline_lty = "dashed",  # 垂直线的线型
                   vertline_col = "black",  # 垂直线的颜色
                   ci_col = "black",           
                   ci_lwd = 2.5,    
                   ci_alpha = 1,
                   point_pch = 22,               
                   point_col = "#E41A1C", 
                   point_fill = "#E41A1C",       
                   point_cex = 1.2,
                   line_gp = gpar(cex = 1.5), 
                   xaxis_lwd = 2)    

pdf("./graph/forest_plot_g130.pdf", width = 6, height = 6)


p <- forest(dt_new[, c("Age", " ", "Interaction")],
            est = dt_new$beta,
            lower = dt_new$CI_lower, 
            upper = dt_new$CI_upper,
            sizes = 0.9,                       
            ci_column = 2,
            vgap = 0.8, 
            ref_line = 0,
            xlim = c(-2.5, 1.),
            ticks_at = c(-2, -1, 0, 1),
            theme = tm,
            footnote = "* p < 0.05; ** p < 0.01",
            # arrow_lab = c("Negative Correlation", "Positive Correlatison"),
            mar = c(4, 4, 2, 2))


p <- edit_plot(p, row = c(gh_row:(gh_row+2)), 
               which = "text",
               col = 1,  
               gp = gpar(col = "#ef8632"))


p <- edit_plot(p, row = c(pe_row:(pe_row+2)), 
               which = "text", 
               col = 1,  
               gp = gpar(col = "#aa6960"))


significant_rows <- which(as.numeric(dt_new$Interaction) < 0.05 & dt_new$Interaction != "")
if(length(significant_rows) > 0) {
  p <- edit_plot(p, row = significant_rows, 
                 which = "text",
                 col = 3,  
                 gp = gpar(col = "#E41A1C"))
}

plot(p)


grid.lines(
  x = c(0, 1),
  y = c(0.05, 0.05), 
  gp = gpar(lwd = 2),
  default.units = "npc"
)


dev.off()









