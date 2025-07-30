rm(list = ls())

library(readxl)
library(writexl)
library(haven)
library(tidyverse)
library(data.table)

library(genepi.utils)

library(CMplot)
library(rMVP)
library(ggplot2)


## fig 4a
pe_meta_5 <- fread('../data/METAANALYSIS1.TBL')%>% 
  mutate(CHR = str_split_i(MarkerName, ":", 1) %>% as.numeric(),
         POS = str_split_i(MarkerName, ":", 2) %>% as.numeric(),
         Allele1 = toupper(Allele1),
         Allele2 = toupper(Allele2))

TSBC_lead_pe <- pe_meta_5 %>% filter(`P-value` < 1e-5)

single_pe <- pe_meta_5 %>% select(MarkerName, Chromosome = CHR, Position = POS, PE = `P-value`) 
MVP.Report(single_pe, plot.type="m",LOG10=TRUE,
           chr.labels=paste("Chr",c(1:22),sep=""),
           threshold=c(1e-5),threshold.lty=c(2), cex.axis = 0.6, cex = 0.6,
           threshold.lwd=c(1), threshold.col=c("grey"),
           col = c("#809bce","#95b8d1","#b8e0d4","#d6eadf","#eac4d5"),
           file.type="pdf", memo="GWASmeta_singlepe5_man",dpi=50,file.output=TRUE,
           verbose=TRUE,width=14,height=6)

