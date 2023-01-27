#Figure3AB and s. Figure2A

# ## tables for heatmap ---------------------------------------------------
# heat map for DEGs --------------------------------------------------------

##### degs_in_time_course analyses
deg_3p.up = ft.deg.cpm.up1
for (i in c(2,3,5,6,7)) {
  deg= get(paste0("ft.deg.cpm.up",i))
  deg_3p.up=rbind(deg_3p.up,deg)
}

deg_3p.dn = ft.deg.cpm.dn1
for (i in c(2,3,5,6,7)) {
  deg= get(paste0("ft.deg.cpm.dn",i))
  deg_3p.dn=rbind(deg_3p.dn,deg)
}

deg.3p.id = unique(c(deg_3p.up$zfin_id_symbol,deg_3p.dn$zfin_id_symbol))

for (i in c(1:8)) {
  deg=get(paste0("deg.cpm",i))
  deg = deg%>% filter(zfin_id_symbol %in% deg.3p.id) %>% arrange(desc(zfin_id_symbol))
  md.fc = deg[,c("mean.cpm.pos","mean.cpm.wt")]
  rownames(md.fc)=deg$zfin_id_symbol
  name1 = paste0("heat.deg.cpm",i)
  assign(name1,md.fc)
}


ht.mt.f3 = data.frame(heat.deg.cpm5,heat.deg.cpm1$mean.cpm.pos,
                      heat.deg.cpm6,heat.deg.cpm2$mean.cpm.pos,
                      heat.deg.cpm7,heat.deg.cpm3$mean.cpm.pos)



colnames(ht.mt.f3)=c("6_bPAC+","6_wt","6_bPAC-","13_bPAC+","13_wt","13_bPAC-",
                     "120_bPAC+","120_wt","120_bPAC-")


# pairwise comparisons
library("PerformanceAnalytics")
pdf(file = "./Figures/pairwise.timepoint.correlation.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8) # The height of the plot in inches

chart.Correlation(as.matrix(t(scale(t(ht.mt.f3)))), histogram=TRUE, pch=19)
dev.off()

#colour code for heatmap
library (gplots)
library(RColorBrewer)

breaks=seq(-3,3,0.01)
#mycol <- colorpanel(n=length(breaksList)-1,low='blue', mid = 'white', high='firebrick3')
#mycol <- colorpanel(1000,"blue","white","red")
mycol = colorRampPalette(c("#313695","#4575B4","#74ADD1","#ABD9E9" ,"#E0F3F8","#FFFFBF" ,"#FEE090" ,"#FDAE61", "#F46D43","#D73027","#A50026" ))(n=600)

#matrix
st_sig_tx =as.matrix(t(scale(t(ht.mt.f3))))
rownames(st_sig_tx)=rownames(ht.mt.f3)

#clustering
distfunc <- function(x) dist(x, method="euclidean")
hclustfunc <- function(x) hclust(x, method="ward.D2")

hc = hclust(dist(t(st_sig_tx), method="euclidean"), method="ward.D2")
dd = as.dendrogram(hc)

#test for optimal number of cluster
#library(factoextra)
#fviz_nbclust(st_sig_tx, FUN = hcut, method = "silhouette")

#clustering
rhc = hclust(dist(st_sig_tx, method="euclidean"), method="ward.D2")
#set number of clusters
gr.row <- cutree(rhc, 4) 
ht.table = data.frame(rhc$labels,rhc$order,gr.row)
#colour for clusters
col1 <- brewer.pal(4, "Set3")
desturate =viridis(4)
#table for clustering
write.csv(ht.table,"./outputs/fig3a.heat.clustering.table.csv")

#clustering
png("./figures/Fig3Atree_heatmap.png", width = 300, height = 200)
plot(dd)
dev.off()

#heatmap
tiff("./figures/Fig3A_3point.heatmap.tiff", width = 12, height = 15, res=1200, units = "cm", compression = "lzw")

heatmap.2((st_sig_tx), col=mycol, scale='none',
               dendrogram = "both",breaks=breaks,hclustfun=hclustfunc, distfun=distfunc,
               #main=c(paste0("Z-score (",nrow(st_sig_tx),")")), 
               Colv = reorder(as.dendrogram(dd),1:9),
               density='none', trace = 'none',
               na.color = F,
               offsetRow = 0,
               offsetCol = c(0,0.2),
               labRow = "",
               colsep = c(4,6,7),
               #rowsep = c(1649,2792,3935), 
               #sepwidth = 2,        
               cexCol = 1,
               #labCol = c("13_pos","6_pos","6_wt","13_wt","LD_wt","120_wt","120_pos","LD_pos"),
               #cellnote=ht.FDR,
               #notecex=1,
               #notecol="black",
               ColSideColors = c(rep(c(desturate[1],desturate[4],desturate[1]),3)),
               RowSideColors=col1[gr.row],
               margins = c(5, 5),
               key.title = "",
               key.xlab = "",
               key.ylab = "",
               keysize = 0.5,
               key = F)
dev.off()

length(rownames(st_sig_tx))

# GO_cluster -------------------------------------------------------
tc_point <- list()
cluster.name = c("cluster1","cluster2","cluster3","cluster4")
for (i in c(1:4)){
  cluster_psites <- which(ht.table$gr.row == i)
  cluster_prots <- unique(ht.table$rhc.labels[cluster_psites])
  tc_point[[cluster.name[i]]] <- cluster_prots
}

#run gprofiler2 externally

#cluster1 =



# venn diagram for DEGs ---------------------------------------------------

###star.pos
up.list = list(
  "up.d6" = ft.deg.cpm.up5$zfin_id_symbol,
  "up.d13" = ft.deg.cpm.up6$zfin_id_symbol,
  "up.d120" = ft.deg.cpm.up7$zfin_id_symbol) 

dn.list = list(
  "down.d6" = ft.deg.cpm.dn5$zfin_id_symbol,
  "down.d13" = ft.deg.cpm.dn6$zfin_id_symbol,
  "down.d120" = ft.deg.cpm.dn7$zfin_id_symbol) 

all.list = list(
  "up.d6" = ft.deg.cpm.up5$zfin_id_symbol,
  "up.d13" = ft.deg.cpm.up6$zfin_id_symbol,
  "up.d120" = ft.deg.cpm.up7$zfin_id_symbol,
  "down.d6" = ft.deg.cpm.dn5$zfin_id_symbol,
  "down.d13" = ft.deg.cpm.dn6$zfin_id_symbol,
  "down.d120" = ft.deg.cpm.dn7$zfin_id_symbol
)


library(eulerr)
library(UpSetR)


fit3= euler(all.list, shape="ellipse")
fit4= euler(up.list, shape="ellipse")
fit5= euler(dn.list, shape="ellipse")


venncol6 = c("#A50026","#F46D43","#FEE090" ,"#313695","#4575B4","#ABD9E9")

deg.all.venn =plot(fit3, 
                   quantities = TRUE,
                   fill = venncol6, alpha =0.3,
                   lty = 1,
                   labels = list(font = 4))

deg.up.venn =plot(fit4, 
                   quantities = TRUE,
                   fill = plasma(6)[4:6], alpha =0.3,
                   lty = 1,
                   labels = list(font = 4))
deg.dn.venn =plot(fit5, 
                   quantities = TRUE,
                   fill = venncol6[4:6], alpha =0.3,
                   lty = 1,
                   labels = list(font = 4))

tiff("./figures/Fig3B_deg.all.venn3.15x12.tiff", width =15 , height = 12, res = 1200, units = "cm",compression = "lzw")
#print(deg.all.venn)
plot_grid(deg.up.venn,deg.dn.venn)
dev.off()
###
'
original fitted residuals regionError
up.d6                                                406    406         0       0.002
up.d13                                               550    550         0       0.003
up.d120                                              527    527         0       0.003
down.d6                                              404    404         0       0.002
down.d13                                             837    837         0       0.004
down.d120                                            508    508         0       0.003
up.d6&up.d13                                         112    112         0       0.001
up.d6&up.d120                                         30     30         0       0.000
up.d6&down.d6                                          0      0         0       0.000
up.d6&down.d13                                        10     10         0       0.000
up.d6&down.d120                                        1      1         0       0.000
up.d13&up.d120                                        69     69         0       0.000
up.d13&down.d6                                         6      0         6       0.001
up.d13&down.d13                                        0      0         0       0.000
up.d13&down.d120                                       1      0         1       0.000
up.d120&down.d6                                        0      0         0       0.000
up.d120&down.d13                                      11      0        11       0.002
up.d120&down.d120                                      0      0         0       0.000
down.d6&down.d13                                     357    357         0       0.002
down.d6&down.d120                                    217    217         0       0.001
down.d13&down.d120                                    84      0        84       0.019
up.d6&up.d13&up.d120                                  48     48         0       0.000
up.d6&up.d13&down.d6                                   0      0         0       0.000
up.d6&up.d13&down.d13                                  0      0         0       0.000
up.d6&up.d13&down.d120                                 0      0         0       0.000
up.d6&up.d120&down.d6                                  0      0         0       0.000
up.d6&up.d120&down.d13                                 1      0         1       0.000
up.d6&up.d120&down.d120                                0      0         0       0.000
up.d6&down.d6&down.d13                                 0      0         0       0.000
up.d6&down.d6&down.d120                                0      0         0       0.000
up.d6&down.d13&down.d120                               0      0         0       0.000
up.d13&up.d120&down.d6                                 0      0         0       0.000
up.d13&up.d120&down.d13                                0      0         0       0.000
up.d13&up.d120&down.d120                               0      0         0       0.000
up.d13&down.d6&down.d13                                0      0         0       0.000
up.d13&down.d6&down.d120                               0      0         0       0.000
up.d13&down.d13&down.d120                              0      0         0       0.000
up.d120&down.d6&down.d13                               1      0         1       0.000
up.d120&down.d6&down.d120                              0      0         0       0.000
up.d120&down.d13&down.d120                             0      0         0       0.000
down.d6&down.d13&down.d120                           277    277         0       0.001
up.d6&up.d13&up.d120&down.d6                           0      0         0       0.000
up.d6&up.d13&up.d120&down.d13                          0      0         0       0.000
up.d6&up.d13&up.d120&down.d120                         0      0         0       0.000
up.d6&up.d13&down.d6&down.d13                          0      0         0       0.000
up.d6&up.d13&down.d6&down.d120                         0      0         0       0.000
up.d6&up.d13&down.d13&down.d120                        0      0         0       0.000
up.d6&up.d120&down.d6&down.d13                         0      0         0       0.000
up.d6&up.d120&down.d6&down.d120                        0      0         0       0.000
up.d6&up.d120&down.d13&down.d120                       0      0         0       0.000
up.d6&down.d6&down.d13&down.d120                       0      0         0       0.000
up.d13&up.d120&down.d6&down.d13                        0      0         0       0.000
up.d13&up.d120&down.d6&down.d120                       0      0         0       0.000
up.d13&up.d120&down.d13&down.d120                      0      0         0       0.000
up.d13&down.d6&down.d13&down.d120                      0      0         0       0.000
up.d120&down.d6&down.d13&down.d120                     0      0         0       0.000
up.d6&up.d13&up.d120&down.d6&down.d13                  0      0         0       0.000
up.d6&up.d13&up.d120&down.d6&down.d120                 0      0         0       0.000
up.d6&up.d13&up.d120&down.d13&down.d120                0      0         0       0.000
up.d6&up.d13&down.d6&down.d13&down.d120                0      0         0       0.000
up.d6&up.d120&down.d6&down.d13&down.d120               0      0         0       0.000
up.d13&up.d120&down.d6&down.d13&down.d120              0      0         0       0.000
up.d6&up.d13&up.d120&down.d6&down.d13&down.d120        0      0         0       0.000

diagError: 0.019 
stress:    0.003 
'

##deg.lis. comparison sp.sn
deg.sp.sn.up.list = list(
  "bPAC-:d6" = c(ft.deg.cpm.up1$zfin_id_symbol),
  "bPAC+:d6" = c(ft.deg.cpm.up5$zfin_id_symbol),
  "bPAC-:d13" = c(ft.deg.cpm.up2$zfin_id_symbol),
  "bPAC+:d13" = c(ft.deg.cpm.up6$zfin_id_symbol),
  "bPAC-:d120" = c(ft.deg.cpm.up3$zfin_id_symbol),
  "bPAC+:d120" = c(ft.deg.cpm.up7$zfin_id_symbol)
)
deg.sp.sn.dn.list = list(
  "bPAC-:d6" = c(ft.deg.cpm.dn1$zfin_id_symbol),
  "bPAC+:d6" = c(ft.deg.cpm.dn5$zfin_id_symbol),
  "bPAC-:d13" = c(ft.deg.cpm.dn2$zfin_id_symbol),
  "bPAC+:d13" = c(ft.deg.cpm.dn6$zfin_id_symbol),
  "bPAC-:d120" = c(ft.deg.cpm.dn3$zfin_id_symbol),
  "bPAC+:d120" = c(ft.deg.cpm.dn7$zfin_id_symbol)
)



library(venn)


tiff("./figures/s.Fig2A_sp.sn.deg.up.venn6.15x12.tiff", width =12 , height = 12, res = 1200, units = "cm",compression = "lzw")
#print(deg.all.venn)
  venn::venn(deg.sp.sn.up.list, ilabels = TRUE, ellipse=T, box=F, zcolor = magma(6), size = 25, cexil = 1.2, cexsn = 1.3)
dev.off()



tiff("./figures/s.Fig2A_sp.sn.deg.dn.venn6.15x12.tiff", width =12 , height = 12, res = 1200, units = "cm",compression = "lzw")
#print(deg.all.venn)
  venn::venn(deg.sp.sn.dn.list, ilabels = TRUE,ellipse=T, box=F, zcolor = viridis(6), size = 25, cexil = 1.2, cexsn = 1.3)
dev.off()



'
                                                           original   fitted residuals regionError
bPAC-:d6                                                         103    0.000   103.000       0.022
bPAC+:d6                                                         594  588.308     5.692       0.010
bPAC-:d13                                                         95    0.000    95.000       0.020
bPAC+:d13                                                        798  794.937     3.063       0.014
bPAC-:d120                                                        18    0.000    18.000       0.004
bPAC+:d120                                                      1006 1002.522     3.478       0.018
bPAC-:d6&bPAC+:d6                                                196  191.971     4.029       0.003
bPAC-:d6&bPAC-:d13                                                12    0.000    12.000       0.003
bPAC-:d6&bPAC+:d13                                                 8    0.000     8.000       0.002
bPAC-:d6&bPAC-:d120                                                0    0.000     0.000       0.000
bPAC-:d6&bPAC+:d120                                                1    0.000     1.000       0.000
bPAC+:d6&bPAC-:d13                                                 8    0.000     8.000       0.002
bPAC+:d6&bPAC+:d13                                               260  278.219   -18.219       0.009
bPAC+:d6&bPAC-:d120                                                3    0.000     3.000       0.001
bPAC+:d6&bPAC+:d120                                              223  244.597   -21.597       0.009
bPAC-:d13&bPAC+:d13                                              541  539.251     1.749       0.010
bPAC-:d13&bPAC-:d120                                               0    0.000     0.000       0.000
bPAC-:d13&bPAC+:d120                                               1    0.000     1.000       0.000
bPAC+:d13&bPAC-:d120                                               0    0.000     0.000       0.000
bPAC+:d13&bPAC+:d120                                             105  135.990   -30.990       0.009
bPAC-:d120&bPAC+:d120                                             27   23.522     3.478       0.000
bPAC-:d6&bPAC+:d6&bPAC-:d13                                        8    0.000     8.000       0.002
bPAC-:d6&bPAC+:d6&bPAC+:d13                                       40   51.073   -11.073       0.003
bPAC-:d6&bPAC+:d6&bPAC-:d120                                       1    0.000     1.000       0.000
bPAC-:d6&bPAC+:d6&bPAC+:d120                                      14   22.032    -8.032       0.002
bPAC-:d6&bPAC-:d13&bPAC+:d13                                      40    0.000    40.000       0.009
bPAC-:d6&bPAC-:d13&bPAC-:d120                                      0    0.000     0.000       0.000
bPAC-:d6&bPAC-:d13&bPAC+:d120                                      0    0.000     0.000       0.000
bPAC-:d6&bPAC+:d13&bPAC-:d120                                      0    0.000     0.000       0.000
bPAC-:d6&bPAC+:d13&bPAC+:d120                                      0    0.000     0.000       0.000
bPAC-:d6&bPAC-:d120&bPAC+:d120                                     0    0.000     0.000       0.000
bPAC+:d6&bPAC-:d13&bPAC+:d13                                      92   99.573    -7.573       0.003
bPAC+:d6&bPAC-:d13&bPAC-:d120                                      0    0.000     0.000       0.000
bPAC+:d6&bPAC-:d13&bPAC+:d120                                      0    0.000     0.000       0.000
bPAC+:d6&bPAC+:d13&bPAC-:d120                                      0    0.000     0.000       0.000
bPAC+:d6&bPAC+:d13&bPAC+:d120                                    224  191.450    32.550       0.003
bPAC+:d6&bPAC-:d120&bPAC+:d120                                    10   31.597   -21.597       0.005
bPAC-:d13&bPAC+:d13&bPAC-:d120                                     0    0.000     0.000       0.000
bPAC-:d13&bPAC+:d13&bPAC+:d120                                    43   58.836   -15.836       0.004
bPAC-:d13&bPAC-:d120&bPAC+:d120                                    0    0.000     0.000       0.000
bPAC+:d13&bPAC-:d120&bPAC+:d120                                    5    0.000     5.000       0.001
bPAC-:d6&bPAC+:d6&bPAC-:d13&bPAC+:d13                             91    0.000    91.000       0.019
bPAC-:d6&bPAC+:d6&bPAC-:d13&bPAC-:d120                             0    0.000     0.000       0.000
bPAC-:d6&bPAC+:d6&bPAC-:d13&bPAC+:d120                             1    0.000     1.000       0.000
bPAC-:d6&bPAC+:d6&bPAC+:d13&bPAC-:d120                             0    0.000     0.000       0.000
bPAC-:d6&bPAC+:d6&bPAC+:d13&bPAC+:d120                            12    9.090     2.910       0.000
bPAC-:d6&bPAC+:d6&bPAC-:d120&bPAC+:d120                            0    0.000     0.000       0.000
bPAC-:d6&bPAC-:d13&bPAC+:d13&bPAC-:d120                            0    0.000     0.000       0.000
bPAC-:d6&bPAC-:d13&bPAC+:d13&bPAC+:d120                            7    0.000     7.000       0.001
bPAC-:d6&bPAC-:d13&bPAC-:d120&bPAC+:d120                           0    0.000     0.000       0.000
bPAC-:d6&bPAC+:d13&bPAC-:d120&bPAC+:d120                           1    0.000     1.000       0.000
bPAC+:d6&bPAC-:d13&bPAC+:d13&bPAC-:d120                            0    0.000     0.000       0.000
bPAC+:d6&bPAC-:d13&bPAC+:d13&bPAC+:d120                           50   41.585     8.415       0.001
bPAC+:d6&bPAC-:d13&bPAC-:d120&bPAC+:d120                           0    0.000     0.000       0.000
bPAC+:d6&bPAC+:d13&bPAC-:d120&bPAC+:d120                           2    0.000     2.000       0.000
bPAC-:d13&bPAC+:d13&bPAC-:d120&bPAC+:d120                          4    0.000     4.000       0.001
bPAC-:d6&bPAC+:d6&bPAC-:d13&bPAC+:d13&bPAC-:d120                   2    0.000     2.000       0.000
bPAC-:d6&bPAC+:d6&bPAC-:d13&bPAC+:d13&bPAC+:d120                  37    0.000    37.000       0.008
bPAC-:d6&bPAC+:d6&bPAC-:d13&bPAC-:d120&bPAC+:d120                  0    0.000     0.000       0.000
bPAC-:d6&bPAC+:d6&bPAC+:d13&bPAC-:d120&bPAC+:d120                  1    0.000     1.000       0.000
bPAC-:d6&bPAC-:d13&bPAC+:d13&bPAC-:d120&bPAC+:d120                 0    0.000     0.000       0.000
bPAC+:d6&bPAC-:d13&bPAC+:d13&bPAC-:d120&bPAC+:d120                 1    0.000     1.000       0.000
bPAC-:d6&bPAC+:d6&bPAC-:d13&bPAC+:d13&bPAC-:d120&bPAC+:d120        0    0.000     0.000       0.000

diagError: 0.022 
stress:    0.014
'


# identification of intersection between bPAC+ and - ----------------------

###list of intersection


#up-regulated DEGs
for (i in c(1:8)) {
  up.l = get(paste0("ft.deg.cpm.up",i))
  up.name = paste0("up.venn",i)
  assign(up.name,up.l$zfin_id_symbol)
}

#down-regulated DEGs
for (i in c(1:8)) {
  dn.l = get(paste0("ft.deg.cpm.dn",i))
  dn.name = paste0("dn.venn",i)
  assign(dn.name,dn.l$zfin_id_symbol)
}

#making lists
up.venn.list = c()
for (i in c(1:8)) {
  up.name = paste0("up.venn",i)
  up.venn.list=c(up.venn.list,up.name)
}

dn.venn.list = c()
for (i in c(1:8)) {
  dn.name = paste0("dn.venn",i)
  dn.venn.list=c(dn.venn.list,dn.name)
}

#up.lists, sn = bPAC-, sp = bPAC+, int = intersection, uq = unique
ven.up.int=list()
ven.up.sn.uq=list()
ven.up.sp.uq=list()
for (i in c(1:4)) {
  sn.up=get(paste0("up.venn",i))
  sp.up=get(paste0("up.venn",i+4))
  upve = list("sn"=sn.up,
              "sp"=sp.up)
  ven.up.int[[i]] = intersect(sn.up,sp.up)
  ven.up.sn.uq[[i]]= sn.up[-c(which(sn.up %in% ven.up.int[[i]]))]
  ven.up.sp.uq[[i]]= sp.up[-c(which(sp.up %in% ven.up.int[[i]]))]
  
}

#down.list,  sn = bPAC-, sp = bPAC+, int = intersection, uq = unique
ven.dn.int=list()
ven.dn.sn.uq=list()
ven.dn.sp.uq=list()
for (i in c(1:4)) {
  sn.dn=get(paste0("dn.venn",i))
  sp.dn=get(paste0("dn.venn",i+4))
  dnve = list("sn"=sn.dn,
              "sp"=sp.dn)
  ven.dn.int[[i]] = intersect(sn.dn,sp.dn)
  ven.dn.sn.uq[[i]]= sn.dn[-c(which(sn.dn %in% ven.dn.int[[i]]))]
  ven.dn.sp.uq[[i]]= sp.dn[-c(which(sp.dn %in% ven.dn.int[[i]]))]
  
}


#intersection

cp.time = c("d6","d13","d120","LD")
int.venn.up = c()
int.venn.dn = c()
tp.up=c()
tp.dn=c()

for (i in 1:4) {
  ven.up = ven.up.int[[i]]
  ven.dn = ven.dn.int[[i]]
  int.venn.up = c(int.venn.up,ven.up)
  int.venn.dn = c(int.venn.dn,ven.dn)
  tp.up=c(tp.up,rep(cp.time[i],length(ven.up)))
  tp.dn=c(tp.dn,rep(cp.time[i],length(ven.dn)))
}

int.venn.up.tb = data.frame("DEGs"=int.venn.up,
                            "timepoint"= tp.up,
                            "regulation" = "up")
int.venn.up.tb = data.frame("DEGs"=int.venn.up,
                            "timepoint"= tp.up,
                            "regulation" = "down")
int.venn.tb = rbind(int.venn.up.tb,int.venn.up.tb)

##s.table1b
write.csv(int.venn.tb, "./outputs/s.table1b_DEGs.vendiagram_sn-sp.csv")
