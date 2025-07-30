library(dplyr)
library(openxlsx)
library(circlize)
library(ggplot2)

current_folder <- './Fig5d/'

split_text <- function(text, max_length = 16) {
  if (nchar(text) <= max_length) {
    return(text)
  } else {
    parts <- strsplit(text, "")[[1]]
    lines <- character()
    current_line <- ""
    for (part in parts) {
      if (nchar(paste0(current_line, part)) + 1 > max_length) {
        lines <- c(lines, current_line)
        current_line <- part
      } else {
        current_line <- paste0(current_line, part)
      }
    }
    lines <- c(lines, current_line)
    return(paste(lines, collapse = "\n"))
  }
}

target <- 'PE'         # 'GH'  for the GH figure in Fig.5d 
direction <- 'forward' 
degs <- list('HDP'=-30, 'PE'=160, 'GH'=-90)
start_deg <- degs[target][[1]]
data <- read.xlsx(paste0(current_folder, './', target, 'V1V1_', direction, '.xlsx'))


dat <- data  %>%  select('from'='name', 'to'='matched_compound', 'value'='prop')
sector_name <- c(data$name, data$matched_compound)
sec_source <- c(data$source, rep('Met', length(dat$to)))

dat$from <- as.character(dat$from)
dat$to <- as.character(dat$to)
arrow_color <- case_match(sec_source, '16S'~'#95C6C6', 'ITS'~'#ebdcb2', 'Species'~'#FFB2B9', 'Genus'~'#7EABCA', 'Met'~'#472D7B3a')

gaps <- data.frame(
    name = sector_name,
    source = sec_source
)
gaps <- unique(gaps)

gaps <- gaps  %>%  mutate(
    color=case_match(source, '16S'~'#95C6C6', 'ITS'~'#ebdcb2', 'Species'~'#FFB2B9', 'Genus'~'#7EABCA', 'Met'~'#472D7B3a')
  )

n_16S <- dim(gaps  %>%  filter(source=='16S'))[1]
n_Taxs <- dim(gaps  %>%  filter(source=='Species'))[1]
n_Taxg <- dim(gaps  %>%  filter(source=='Genus'))[1]

bgap <- 5
gaps$gap <- 1
gaps$gap[n_16S] <- bgap       
gaps$gap[n_16S+n_Taxs] <- bgap 
gaps$gap[n_16S+n_Taxs+n_Taxg] <- bgap
gaps$gap[length(gaps$gap)] <- bgap

gap <- setNames(gaps$gap, gaps$name)
grid.col <- setNames(gaps$color, gaps$name)
group <- setNames(gaps$source, gaps$name)
sector_abbr <- data$compound_abbr


pdf(paste0(current_folder, "Fig.5d.pdf"), width = 15, height = 15)

circos.clear()
circos.par(start.degree = start_deg, 
           gap.degree = 2, 
           points.overflow.warning = FALSE
          #  gap.after=gap
           )

chordDiagram(dat,
             grid.col = grid.col,
             col = arrow_color,
             group = group,
             transparency = 0.25,
             directional = 1,
             direction.type = c("arrows", "diffHeight"), 
             diffHeight  = -0.04,
             annotationTrack = "grid", 
             annotationTrackHeight = c(0.08, 0.1),
             link.arr.type = "big.arrow", 
            #  link.sort = TRUE, 
            #  link.decreasing = TRUE,
             link.largest.ontop = TRUE,
             preAllocateTracks = list(
               track.height = 0.1,
               track.margin = c(0.01, 0)
             ))

circos.trackPlotRegion(
  track.index = 2, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    sector.index = sapply(CELL_META$sector.index, split_text)

    circos.text(
      x = CELL_META$xcenter,  
      y = 1.2,                
      labels = sector.index, 
      facing = "clockwise",   
      niceFacing = TRUE, 
      cex = 1,        
      col = "black",  
      adj = c(0, 0.5))   
  }
)


circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    circos.text(
      x = CELL_META$xcenter, 
      y = 3.5,                
      labels = rep(' ', 10), 
      facing = "clockwise",   
      niceFacing = TRUE, 
      cex = 1,        
      col = "black",  
      adj = c(0, 0.5))
  }
)

dev.off()
