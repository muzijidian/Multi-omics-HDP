
rm(list = ls())
setwd('/Users/mzjd/Documents/HDP-multiomics/haonan_code')

library(grid)
library(forestploter)
library(readxl)

dt = read_excel("../Anti_data/G130_forest.xlsx")


dt$Age[dt$Age == "GH (48/18)"] <- "GH (T1)"
dt$Age[dt$Age == "< 30 (1039)"] <- "< 30 (48/1,039)"
dt$Age[dt$Age == "≥ 30 (295)" | dt$Age == "... 30 (295)"] <- "\u2265 30 (18/295)" 

dt$Age[dt$Age == "PE (27/11)"] <- "PE (T2)"
dt$Age[dt$Age == "< 30 (790)"] <- "< 30 (27/790)"
dt$Age[dt$Age == "≥ 30 (210)" | dt$Age == "... 30 (210)"] <- "\u2265 30 (11/210)"


dt$CI_lower <- dt$beta - 1.96 * dt$se
dt$CI_upper <- dt$beta + 1.96 * dt$se


dt$` ` <- paste(rep(" ", 15), collapse = " ")

dt$Interaction <- ifelse(!is.na(dt$interaction), 
                         sprintf("%.3f", dt$interaction), 
                         "")


tm <- forest_theme(base_size = 16,
                   core = list(padding = unit(c(7, 11), "mm")), 
                   refline_gp = gpar(lwd = 2, lty = "dashed", col = "#656c8a"),
                   vertline_lwd = 2, vertline_lty = "dashed", vertline_col = "black",
                   ci_col = "black", ci_lwd = 2.5, ci_alpha = 1,
                   point_pch = 22, point_col = "#E41A1C", point_fill = "#E41A1C", point_cex = 1.2,
                   line_gp = gpar(cex = 1.5), 
                   xaxis_gp = gpar(lwd = 2, cex = 1.5))


quartz(type = "pdf", file = "./Anti_graph/fig3f_forest_plot_g130.pdf", width = 7, height = 6)

gh_row <- which(dt$Age == "GH (T1)")
pe_row <- which(dt$Age == "PE (T2)")

p <- forest(dt[, c("Age", " ", "Interaction")],
            est = dt$beta,
            lower = dt$CI_lower, 
            upper = dt$CI_upper,
            sizes = 0.9,                       
            ci_column = 2,
            vgap = 0.8, 
            ref_line = 0,
            xlim = c(-2.5, 1.),
            ticks_at = c(-2, -1, 0, 1),
            theme = tm,
            mar = c(4, 4, 2, 2))

p <- edit_plot(p, row = c(gh_row:(gh_row+2)), 
               which = "text",
               col = 1,  
               gp = gpar(col = "#ef8632"))

p <- edit_plot(p, row = c(pe_row:(pe_row+2)), 
               which = "text", 
               col = 1,  
               gp = gpar(col = "#aa6960"))


significant_rows <- which(as.numeric(dt$Interaction) < 0.05 & dt$Interaction != "")
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






