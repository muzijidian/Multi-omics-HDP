library(dplyr)
library(tidyr)
library(openxlsx)
library(circlize)
library(ggplot2)
library(stringr)

# split long text
split_text <- function(text, max_length = 32) {
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


Draw_chrod <- function(dat, gaps, save_file, start_deg=85, hw=20){
  gaps <- unique(gaps)
  gap <- setNames(gaps$gap, gaps$name)
  grid.col <- setNames(gaps$color, gaps$name)
  group <- gaps  %>%  select(name, source)  %>%  arrange(source)
  group <- setNames(group$source, group$name)
  
  arrow_color <- dat  %>%  left_join(gaps, by=c('from'='name'))  %>%  select(color)
  arrow_color <- arrow_color[[1]]

  pdf(save_file, width = hw, height = hw)
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
  print(save_file)
}


Prepare_Micros <- function(file, thresh=0.05){
    data <- read.csv(file)
    data <- data  %>%  filter(b1_pv_FDR<thresh)
    data <- data  %>%  arrange(name.x, name.y)

    dat <- data  %>%  select('from'=name.x, 'to'=name.y, 'value'='b1_est')
    dat$from <- as.character(dat$from)
    dat$to <- as.character(dat$to)

    sector_name <- c(data$name.x, data$name.y)
    # sec_source <- c(data$source, data$source)
    sec_source <- c(data$source, rep('ITS', length(data$name.y))) 

    gaps <- data.frame(
        name = sector_name,
        source = sec_source
    )
    gaps <- unique(gaps)

    gaps <- gaps  %>%  mutate(
        source=case_when(source=='16S'~'16S',
                          source=='ITS'~'ITS',
                          source=='genus'~'Tax_g',
                          source=='species'~'Tax_s')
        )

    gaps <- gaps  %>%  mutate(
        color=case_match(source, '16S'~'#95C6C6', 'ITS'~'#ebdcb2', 'Tax_s'~'#FFB2B9', 'Tax_g'~'#7EABCA', 'Met'~'#472D7B3a')
        )

    n_micro_x <- dim(gaps  %>%  filter(source=='Tax_s'))[1]
    n_micro_y <- dim(gaps  %>%  filter(source=='16S'))[1]

    bgap <- 5
    gaps$gap <- 1
    gaps$gap[n_micro_x] <- bgap       
    gaps$gap[n_micro_x+n_micro_y] <- bgap 

    return(list('dat'=dat, 'gaps'=gaps))
}


file <- paste0('./data/sanky_bac_fun_V12_antibiotic_ascov.csv')
DATA <- Prepare_Micros(file, 0.05)
dat <- DATA$dat
gaps <- DATA$gaps

save_file <- paste0('./Fig3d,3e/Fig.3d.pdf')
Draw_chrod(dat, gaps, save_file, start_deg=125)


# file <- paste0('./HDP_CLPM/sanky_bac_V13_antibiotic_ascov.csv')
# DATA <- Prepare_Micros(file, 0.1)
# dat <- DATA$dat
# gaps <- DATA$gaps

# save_file <- paste0('./HDP_CLPM/Fig.3e.pdf')
# Draw_chrod(dat, gaps, save_file)
