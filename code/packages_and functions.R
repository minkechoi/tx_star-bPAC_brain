###package and functions 

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.

# load some packages
library(tidyverse)
library(cowplot)
library(png)
library(hrbrthemes)
library(viridis)
library(limma) 
library(edgeR) 
library(dplyr)
library(genefilter)
library(devtools)
library(gplots)
library(ggplot2)
library(calibrate)
library(graphics)
library(data.table)
library(clusterProfiler)
library(factoextra)
library(stringr)
library(RColorBrewer)
options("scipen" = 100)


#define some colors
#From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
               "#CC79A7", "#000000")
########################## SHOW as multiple PIE Charts ###################
pie(rep(1,8), col=Okabe_Ito, Okabe_Ito, main='Okabe Ito')

#others

show_colors <- function(colors) { 
  ggplot(data.frame(id=seq_along(colors), color=colors)) + 
    geom_tile(aes(id, 1, fill=color)) + 
    scale_fill_identity()+
    theme_minimal()+
    theme(axis.text = element_text(size = 0))+
    labs(x="",y="")
}


#save session info and Rstudio version info for reproducibility
writeLines(capture.output(sessionInfo()), "code/sessionInfo.txt")
writeLines(capture.output(rstudioapi::versionInfo()), "code/versionInfo.txt")
