#helen's
library(tidyverse)
library(factoextra)
library(stringr)
library(viridis)
library(clusterProfiler)

#read tables
for (i in 5:7) {
  files=paste0("./data/topqlf.cpm.",i,".csv")
  deg= read.csv(files,row.names = 1)
  name = paste0("deg.cpm_",i)
  assign(name, deg)
  print(i)
}

time.name = c("d6.sn", "d13.sn", "d120.sn", "LD.sn","d6", "d13", "d120", "LD")

for(i in 5:7){
  tb = get(paste0("deg.cpm_",i))
  ttb=c()
  mcp=c()
  mcw=c()
  colnames(tb)=c("logFC","logCPM","F","PValue","FDR","id","ensembl_gene_id","zfin_id_symbol","description","pos1","pos2","pos3","pos4","pos5","wt1","wt2","wt3","wt4","wt5")
  tb = tb %>% arrange(desc(logCPM))
  d = duplicated(tb$zfin_id_symbol)
  tb=tb[!d,]
  
  for (k in 1:nrow(tb)) {
    tt=t.test(as.matrix(tb[k,c(10:14)]),as.matrix(tb[k,c(15:19)]))
    mcp=c(mcp,tt$estimate[1])
    mcw=c(mcw,tt$estimate[2])
    ttb=c(ttb,tt$p.value)
  }
  tb[,"mean.cpm.pos"]=mcp
  tb[,"mean.cpm.wt"]=mcw
  tb[,"t.test"]=ttb
  tb[,"source"] = time.name[i]
  name = paste0("deg.cpm_",i)
  assign(name,tb)
}

gene_all = deg.cpm_5
for (i in c(6:7)) {
  deg= get(paste0("deg.cpm_",i))
  gene_all=rbind(gene_all,deg)
}


prot_id_map <- clusterProfiler::bitr(unique(gene_all$zfin_id_symbol), "SYMBOL",
                                     "ENTREZID", "org.Dr.eg.db")
d=duplicated(prot_id_map$SYMBOL)
prot_id_map=prot_id_map[!d,]
rownames(prot_id_map) <- prot_id_map$SYMBOL



#list of background genes
for (i in c(5:7)){
  deg = get(paste0("deg.cpm_",i))
  deg = deg %>% dplyr::filter(mean.cpm.pos > 5|mean.cpm.wt >5)
  deg = prot_id_map[deg$zfin_id_symbol, "ENTREZID"]
  deg = na.omit(deg)
  name = paste0("bk.gene",i)
  assign(name, deg)
  print(i)
}

bk.all.basal = unique(c(bk.gene5,bk.gene6,bk.gene7))


###degs

#seperating DEGs


for(i in 5:7){
  deg=get(paste0("deg.cpm_",i))
  deg = deg %>% filter(mean.cpm.pos > 5|mean.cpm.wt >5)%>% filter(t.test < 0.05)%>% 
    filter(FDR < 0.05)%>%filter(zfin_id_symbol %in% prot_id_map$SYMBOL) %>% arrange(desc(logCPM))
  d = duplicated(deg$zfin_id_symbol)
  deg=deg[!d,]
  name = paste0("ft.deg.cpm",i)
  assign(name,deg)
}

#listing of significant DEGs 

b.ft.id = ft.deg.cpm5$id
for (i in c(6:7)) {
  deg= get(paste0("ft.deg.cpm",i))
  b.ft.id=unique(c(b.ft.id,deg$id))
}

#ft_gene_all
ft.gene_all = ft.deg.cpm5
for (i in c(6:7)) {
  deg= get(paste0("ft.deg.cpm",i))
  ft.gene_all=rbind(ft.gene_all,deg)
}

for (i in c(5:7)) {
  deg=get(paste0("ft.deg.cpm",i))
  deg.up = deg %>% filter(logFC > 1) %>% filter(FDR < 0.01)%>% arrange(desc(logFC))
  deg.dn = deg %>% filter(logFC < -1) %>% filter(FDR < 0.01) %>% arrange(logFC)
  name1 = paste0("ft.deg.cpm.up",i)
  assign(name1,deg.up)
  name2 = paste0("ft.deg.cpm.dn",i)
  assign(name2,deg.dn)
  
}

# adjust of FC threshold to 1.5FC for 120 dpf 

for (i in c(7)) {
  deg=get(paste0("ft.deg.cpm",i))
  deg.up = deg %>% filter(logFC > log2(1.5)) %>% filter(FDR < 0.01)%>% arrange(desc(logFC))
  deg.dn = deg %>% filter(logFC < -log2(1.5))  %>% filter(FDR < 0.01)%>% arrange(logFC)
  name1 = paste0("ft.deg.cpm.up",i)
  assign(name1,deg.up)
  name2 = paste0("ft.deg.cpm.dn",i)
  assign(name2,deg.dn)
  
}

# combining tables 
deg_all.up = ft.deg.cpm.up5
for (i in c(6:7)) {
  deg= get(paste0("ft.deg.cpm.up",i))
  deg_all.up=rbind(deg_all.up,deg)
}

deg_all.dn = ft.deg.cpm.dn5
for (i in c(6:7)) {
  deg= get(paste0("ft.deg.cpm.dn",i))
  deg_all.dn=rbind(deg_all.dn,deg)
}

deg.all = rbind(deg_all.up,deg_all.dn)


#down degs
tc_point_dn <- list()
## Reformat the timecourse cluster classifications as a list

for (i in c(5:7)){
  cluster_psites <- which(deg_all.dn$source == time.name[i])
  cluster_prots <- unique(deg_all.dn$zfin_id_symbol[cluster_psites])
  cluster_prots_entrez <- prot_id_map[cluster_prots, "ENTREZID"]
  cluster_prots_entrez=na.omit(cluster_prots_entrez)
  tc_point_dn[[time.name[i]]] <- cluster_prots_entrez
}
## Do the enrichment analysis
##GOBP
tc_enrich.dn <- clusterProfiler::compareCluster(tc_point_dn, fun = "enrichGO",  OrgDb='org.Dr.eg.db', ont = "BP",
                                                bk.all.basal,pvalueCutoff = 0.05,
                                                pAdjustMethod = "BH",qvalueCutoff = 0.2,
                                                minGSSize = 10,
                                                maxGSSize = 10000, readable = TRUE)

##tables

tc_enrich.dn.bp.table = data.frame(tc_enrich.dn[1:length(tc_enrich.dn[]$ID)])

##filter neurogenesis
neurogenesis=dplyr::filter(tc_enrich.dn.bp.table, Description == "neurogenesis")


d6.neurogenesis= str_split(neurogenesis[1,9],pattern = "/")
d6.neurogenesis= data.frame("gene"=Reduce(rbind,d6.neurogenesis))

d13.neurogenesis= str_split(neurogenesis[2,9],pattern = "/")
d13.neurogenesis= data.frame("gene"=Reduce(rbind,d13.neurogenesis))

d120.neurogenesis= str_split(neurogenesis[3,9],pattern = "/")
d120.neurogenesis= data.frame("gene"=Reduce(rbind,d120.neurogenesis))

neurogen = unique(c(d6.neurogenesis$gene,d13.neurogenesis$gene,d120.neurogenesis$gene))

####venn diagram

neurogen.list = list(
  "d6" = d6.neurogenesis$gene,
  "d13" = d13.neurogenesis$gene,
  "d120" = d120.neurogenesis$gene
)


library(eulerr)
library(UpSetR)


fit3= euler(neurogen.list, shape="ellipse")

venncol3 = c("#A50026","#F46D43","#ABD9E9")

deg.neurogen.venn =plot(fit3, 
                   quantities = TRUE,
                   fill = venncol3, alpha =0.3,
                   lty = 1,
                   labels = list(font = 4))


tiff("./figures/neurogenesis.deg.tiff", width =10 , height = 8, res = 1200, units = "cm",compression = "lzw")
#print(deg.all.venn)
plot(deg.neurogen.venn)
dev.off()

####heatmap

#d6,13
neurogen2 = unique(c(d6.neurogenesis$gene,d13.neurogenesis$gene))

#mean_CPM
for (i in c(5:6)) {
  deg=get(paste0("deg.cpm_",i))
  deg = deg%>% filter(zfin_id_symbol %in% neurogen2) %>% arrange(desc(zfin_id_symbol))
  md.fc = deg[,c("mean.cpm.pos","mean.cpm.wt")]
  rownames(md.fc)=deg$zfin_id_symbol
  name1 = paste0("heat.deg.cpm2_",i)
  assign(name1,md.fc)
}

ht.mt.f2 = data.frame(heat.deg.cpm2_5,
                      heat.deg.cpm2_6)

rownames(ht.mt.f2)=deg$zfin_id_symbol
colnames(ht.mt.f2)=c("6_bPAC+","6_wt","13_bPAC+","13_wt")


#colour code for heatmap
library (gplots)
library(RColorBrewer)

breaks=seq(-3,3,0.01)
#mycol <- colorpanel(n=length(breaksList)-1,low='blue', mid = 'white', high='firebrick3')
#mycol <- colorpanel(1000,"blue","white","red")
mycol = colorRampPalette(c("#313695","#4575B4","#74ADD1","#ABD9E9" ,"#E0F3F8","#FFFFBF" ,"#FEE090" ,"#FDAE61", "#F46D43","#D73027","#A50026" ))(n=600)

#matrix
st_sig_tx =as.matrix(t(scale(t(ht.mt.f2))))
#st_sig_tx =as.matrix(ht.mt.f2)

rownames(st_sig_tx)=rownames(ht.mt.f2)

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
gr.row <- cutree(rhc, 6) 
ht.table = data.frame(rhc$labels,rhc$order,gr.row)
#colour for clusters
col1 <- brewer.pal(6, "Set3")
desturate =viridis(4)
#table for clustering
write.csv(ht.table,"./outputs/neurogenesis_heatmap.table.csv")

#clustering
#png("./figures/Fig3Atree_heatmap.png", width = 300, height = 200)
plot(dd)
#dev.off()

#heatmap
tiff("./figures/neurogenesis_heatmap_d613.tiff", width = 6, height = 10, res=1200, units = "cm", compression = "lzw")

heatmap.2((st_sig_tx), col=mycol, scale='none',
          dendrogram = "both",breaks=breaks,hclustfun=hclustfunc, distfun=distfunc,
          #main=c(paste0("Z-score (",nrow(st_sig_tx),")")), 
          Colv = reorder(as.dendrogram(dd),1:9),
          density='none', trace = 'none',
          na.color = F,
          offsetRow = 0,
          offsetCol = c(0,0.2),
          labRow = "",
          colsep = c(2),
          #rowsep = c(1649,2792,3935), 
          #sepwidth = 2,        
          cexCol = 1,
          #labCol = c("13_pos","6_pos","6_wt","13_wt","LD_wt","120_wt","120_pos","LD_pos"),
          #cellnote=ht.FDR,
          #notecex=1,
          #notecol="black",
          #ColSideColors = c(rep(c(desturate[1],desturate[2]),2)),
          RowSideColors=col1[gr.row],
          margins = c(5, 5),
          key.title = "",
          key.xlab = "",
          key.ylab = "",
          keysize = 0.5,
          key = F)
dev.off()

length(rownames(st_sig_tx))

#FC

for (i in c(5:7)) {
  deg=get(paste0("deg.cpm_",i))
  deg = deg%>% filter(zfin_id_symbol %in% neurogen2) %>% arrange(desc(zfin_id_symbol))
  md.fc = deg[,"logFC"]
  name1 = paste0("heat.deg.cpm2_",i)
  assign(name1,md.fc)
}

ht.mt.f2 = data.frame(heat.deg.cpm2_5,
                      heat.deg.cpm2_6)

ht.mt.f3 = data.frame(heat.deg.cpm2_5,
                      heat.deg.cpm2_6,
                      heat.deg.cpm2_7)

rownames(ht.mt.f2)=deg$zfin_id_symbol
rownames(ht.mt.f3)=deg$zfin_id_symbol


colnames(ht.mt.f2)=c("6_bPAC+","13_bPAC+")
colnames(ht.mt.f3)=c("6_bPAC+","13_bPAC+","120_bPAC+")


#colour code for heatmap
library (gplots)
library(RColorBrewer)

breaks=seq(-3,3,0.01)
#mycol <- colorpanel(n=length(breaksList)-1,low='blue', mid = 'white', high='firebrick3')
#mycol <- colorpanel(1000,"blue","white","red")
mycol = colorRampPalette(c("#313695","#4575B4","#74ADD1","#ABD9E9" ,"#E0F3F8","#FFFFBF" ,"#FEE090" ,"#FDAE61", "#F46D43","#D73027","#A50026" ))(n=600)

#matrix
#st_sig_tx =as.matrix(t(scale(t(ht.mt.f2))))
st_sig_tx =as.matrix(ht.mt.f2)

rownames(st_sig_tx)=rownames(ht.mt.f2)
  #for 3 time points
  st_sig_tx =as.matrix(ht.mt.f3)

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
gr.row <- cutree(rhc, 6)
#manual matching between cluster and heatmap
gr.row[gr.row == 1] <-"a"
gr.row[gr.row == 2] <-"b"
gr.row[gr.row == 3] <-"c"
gr.row[gr.row == 4] <-"d"
gr.row[gr.row == 5] <-"e"
gr.row[gr.row == 6] <-"f"

gr.row[gr.row == "a"] <-4
gr.row[gr.row == "b"] <-2
gr.row[gr.row == "c"] <-3
gr.row[gr.row == "d"] <-6
gr.row[gr.row == "e"] <-1
gr.row[gr.row == "f"] <-5  
  
ht.table = data.frame(rhc$labels,rhc$order,gr.row)

#colour for clusters
col1 <- brewer.pal(6, "Set3")
desturate =viridis(2)
#table for clustering
write.csv(ht.table,"./outputs/neurogenesis_heatmap_FC.table.csv")

#clustering
#png("./figures/Fig3Atree_heatmap.png", width = 300, height = 200)
plot(dd)
#dev.off()

tiff("./figures/neurogenesis_heatmap_d613_FC.tiff", width = 7, height = 15, res=1200, units = "cm", compression = "lzw")

heatmap.2((st_sig_tx), col=mycol, scale='none',
          dendrogram = "row",breaks=breaks,hclustfun=hclustfunc, distfun=distfunc,
          #main=c(paste0("Z-score (",nrow(st_sig_tx),")")), 
          #Colv = reorder(as.dendrogram(dd),1:9),
          #Rowv = reorder(as.dendrogram(dd),1:225),
          density='none', trace = 'none',
          na.color = F,
          offsetRow = 0,
          offsetCol = c(0,0.2),
          #labRow = "",
          #colsep = c(1),
          #rowsep = c(1649,2792,3935), 
          #sepwidth = 2,        
          cexCol = 2,
          labCol = c("6_pos","13_pos","120_pos"),
          #cellnote=ht.FDR,
          #notecex=1,
          #notecol="black",
          #ColSideColors = c(rep(c(desturate[1],desturate[2]),2)),
          RowSideColors=col1[as.numeric(gr.row)],
          margins = c(8, 8),
          key.title = "",
          key.xlab = "",
          key.ylab = "",
          keysize = 0.5,
          key = F,
          lwid=c(1,2))
dev.off()

# 3points, heat:rhc= 1,2,3,4,5,6,:2,3,1,5,6,4 
tiff("./figures/neurogenesis_heatmap_d613120_FC.tiff", width = 8, height = 15, res=1200, units = "cm", compression = "lzw")

heatmap.2((st_sig_tx), col=mycol, scale='none',
          dendrogram = "row",breaks=breaks,hclustfun=hclustfunc, distfun=distfunc,
          #main=c(paste0("Z-score (",nrow(st_sig_tx),")")), 
          #Colv = reorder(as.dendrogram(dd),1:9),
          #Rowv = reorder(as.dendrogram(dd),1:225),
          density='none', trace = 'none',
          na.color = F,
          offsetRow = 0,
          offsetCol = c(0,0.2),
          labRow = "",
          #colsep = c(1),
          #rowsep = c(1649,2792,3935), 
          #sepwidth = 2,        
          cexCol = 2,
          labCol = c("6_pos","13_pos","120_pos"),
          #cellnote=ht.FDR,
          #notecex=1,
          #notecol="black",
          #ColSideColors = c(rep(c(desturate[1],desturate[2]),2)),
          RowSideColors=col1[as.numeric(gr.row)],
          margins = c(8, 8),
          key.title = "",
          key.xlab = "",
          key.ylab = "",
          keysize = 0.5,
          key = F,
          lwid=c(1,2))
dev.off()

length(rownames(st_sig_tx))

#cluster selection
ng_c1=ht.table$rhc.labels[c(which(ht.table$gr.row == 1))]
ng_c2=ht.table$rhc.labels[c(which(ht.table$gr.row == 2))]
ng_c3=ht.table$rhc.labels[c(which(ht.table$gr.row == 3))]
ng_c4=ht.table$rhc.labels[c(which(ht.table$gr.row == 4))]
ng_c5=ht.table$rhc.labels[c(which(ht.table$gr.row == 5))]
ng_c6=ht.table$rhc.labels[c(which(ht.table$gr.row == 6))]

#Welch Two Sample t-test
d6fcs=c()
d13fcs=c()
ttb=c()
tt.ste=c()
sd.6=c()
sd.13=c()
for (i in 1:6) {
  deg=ht.mt.f2
  cl=get(paste0("ng_c",i))
  d6fc=deg[c(which(row.names(deg)%in%cl)),1]
  d13fc=deg[c(which(row.names(deg)%in%cl)),2]
  tt=t.test(as.matrix(d6fc),as.matrix(d13fc))
  d6fcs=c(d6fcs,tt$estimate[1])
  d13fcs=c(d13fcs,tt$estimate[2])
  ttb=c(ttb,tt$p.value)
  tt.ste=c(tt.ste,tt$stderr)
  sd.6=c(sd.6,sd(d6fc))
  sd.13=c(sd.13,sd(d13fc))
  deg["cluster"]=paste0("cluster_",i)
  assign(paste0("neurog.deg.c",i),deg)

}

neurog.deg=rbind(neurog.deg.c1,neurog.deg.c2,neurog.deg.c3,neurog.deg.c4,neurog.deg.c5,neurog.deg.c6)

fc_tt = data.frame("d6_mean"=d6fcs,
                     "d13_mean"=d13fcs,
                     "p.val"=ttb,
                     'ste'=tt.ste,
                    "sd.6"=sd.6,
                   "sd.13"=sd.13
                   )

row.names(fc_tt)=c("cluster 1","cluster 2","cluster 3","cluster 4","cluster 5","cluster 6")


fc_tt_long <- data.frame("mean_logFC"=c(fc_tt$d6_mean,fc_tt$d13_mean),
                         "time"= c(rep("6",6),rep("13",6)),
                         "cluster"= rep(c("cluster 1","cluster 2","cluster 3","cluster 4","cluster 5","cluster 6"),2),
                         "sd"=c(fc_tt$sd.6,fc_tt$sd.13)
                         )

fc_tt_long$time <- factor(fc_tt_long$time, levels = c("6","13")) 

tiff("./figures/neurogenesis_d613_mean_FC.tiff", width = 12, height = 15, res=1200, units = "cm", compression = "lzw")

ggplot(data=fc_tt_long, aes(x=time, y=mean_logFC, group=cluster)) +
  geom_line(aes(color=cluster), size=1.5,position=position_dodge(0.2))+
  geom_point(aes(fill=cluster),color="black",shape=22,size=3,position=position_dodge(0.2))+
  geom_errorbar(aes(ymin=mean_logFC-sd, ymax=mean_logFC+sd, color=cluster), width=.2,position=position_dodge(0.2)) +
  theme_classic()+ scale_y_continuous(limits =c(-4.5,0.5))+ 
  scale_color_manual(values = col1)+
  scale_fill_manual(values = col1)+
  geom_hline(yintercept=-log2(1.5), linetype="dashed", color = "grey", size=0.8)+
  #scale_linetype_manual(values=c("dotted", "dotted","dotted","twodash","dotted","twodash"))+
  theme(legend.title = element_blank()) 
dev.off()



for (k in 1:nrow(tb)) {
  tt=t.test(as.matrix(tb[k,c(10:14)]),as.matrix(tb[k,c(15:19)]))
  mcp=c(mcp,tt$estimate[1])
  mcw=c(mcw,tt$estimate[2])
  ttb=c(ttb,tt$p.value)
}
tb[,"mean.cpm.pos"]=mcp
tb[,"mean.cpm.wt"]=mcw
tb[,"t.test"]=ttb

####trends plot for cluster 4 and 6


###selet cpm
#mean_CPM
for (i in c(5:6)) {
  deg=get(paste0("deg.cpm_",i))
  deg = deg%>% filter(zfin_id_symbol %in% neurogen2) %>% arrange(desc(zfin_id_symbol))
  md.fc = deg[,c("mean.cpm.pos","mean.cpm.wt")]
  rownames(md.fc)=deg$zfin_id_symbol
  name1 = paste0("heat.deg.cpm2_",i)
  assign(name1,md.fc)
}

ht.mt.f2 = data.frame(heat.deg.cpm2_5,
                      heat.deg.cpm2_6)

rownames(ht.mt.f2)=deg$zfin_id_symbol
colnames(ht.mt.f2)=c("d6_pos_cpm","d6_wt_cpm","d13_pos_cpm","d13_wt_cpm")
cpm.data.scaled= t(scale(t(as.matrix(ht.mt.f2)),center = F,scale = T))
wt.cpm.data.scaled =data.frame("gene"=rownames(ht.mt.f2),cpm.data.scaled[,c(2,4)],"cluster"=ht.table$gr.row)
pos.cpm.data.scaled =data.frame("gene"=rownames(ht.mt.f2),cpm.data.scaled[,c(1,3)],"cluster"=ht.table$gr.row)

time.name =c("d6_wt_cpm","d13_wt_cpm")
time.name2 = c("d6_pos_cpm","d13_pos_cpm")





## We only need the quantitative data for clustering and it should be
## in a matrix format, not a data.frame (or tibble)

tc_data_long <- wt.cpm.data.scaled %>%
  tidyr::pivot_longer(cols = d6_wt_cpm:d13_wt_cpm, names_to = "sample", values_to = "cpm") %>%
  tidyr::unite(trace, gene, sep="-", remove = FALSE) %>%
  group_by(trace)

tc_data_long2 <- pos.cpm.data.scaled %>%
  tidyr::pivot_longer(cols = d6_pos_cpm:d13_pos_cpm, names_to = "sample", values_to = "cpm") %>%
  tidyr::unite(trace, gene, sep="-", remove = FALSE) %>%
  group_by(trace)


tc_data_long$sample <- factor(tc_data_long$sample, levels = time.name) 
#timepoint
tim= c(1,2)
for(i in 1:2){
  tc_data_long[c(which(tc_data_long$sample == time.name[i] )),"timepoint"]= tim[i]
}

tc_data_long2$sample <- factor(tc_data_long2$sample, levels = time.name2) 
#timepoint
tim= c(1,2)
for(i in 1:2){
  tc_data_long2[c(which(tc_data_long2$sample == time.name2[i] )),"timepoint"]= tim[i]
}


## Make the cluster labels more readable for a figure
tc_data_long$cluster <- factor(tc_data_long$cluster,
                               levels = sort(unique(tc_data_long$cluster)),
                               labels = paste("cluster", sort(unique(tc_data_long$cluster))))

tc_data_long2$cluster <- factor(tc_data_long2$cluster,
                                levels = sort(unique(tc_data_long2$cluster)),
                                labels = paste("cluster", sort(unique(tc_data_long2$cluster))))


#write.csv(tc_data_long,"./TC/tc.mcluster30.5points.csv")

## Make an auxiliary table with mean log2 fold changes for each
## cluster.  We can add these to the plot for easier visualisation.
tc_data_sum <- dplyr::group_by(tc_data_long, cluster, timepoint) %>%
  dplyr::summarise(mean.cpm = mean(cpm))
tc_data_sum2 <- dplyr::group_by(tc_data_long2, cluster, timepoint) %>%
  dplyr::summarise(mean.cpm = mean(cpm))
#write.csv(tc_data_sum,"./TC/tc.mcluster9.3points.sum.csv")
## Produce a plot with one panel per cluster, showing a trace for each
## timecourse for each stimulus.  Also include a loess-smoothed mean
## trace.



tc_data_long.wt=tc_data_long
tc_data_long.pos=tc_data_long2
tc_data_long.wt[,"genotype"]="wt"
tc_data_long.pos[,"genotype"]="star:bPAC"

tc_data_long.wt.pos = rbind(tc_data_long.wt,tc_data_long.pos)



write.csv(tc_data_long.wt.pos,"./outputs/tc_data_long12.wt.pos.csv")
#all
tiff(file = "./figures/neurogenesis_all.tiff",   # The directory you want to save the file in
     width = 10, # The width of the plot in inches
     height = 10, res=1200, units = "cm", compression = "lzw")
# The height of the plot in inches

#tc_data_long.wt.pos$cluster <- factor(tc_data_long.wt.pos$cluster,      # Reordering group factor levels
#                                      levels = c("cluster 2","cluster 7","cluster 4","cluster 8","cluster 1","cluster 3","cluster 5","cluster 6","cluster 9"))
tc_clust_plot <- ggplot2::ggplot(tc_data_long.wt.pos, ggplot2::aes(timepoint, cpm, colour = genotype)) +
  ggplot2::geom_hline(color = "gray", yintercept = 0, linetype = 2, size = 0.5) +
  #ggplot2::geom_line(ggplot2::aes(group = trace), size = 1, alpha = 0.05) +
  ggplot2::geom_boxplot(ggplot2::aes(group = sample), size = 0.5,alpha = 0.1,outlier.size = 0.5) +
  ggplot2::geom_point(ggplot2::aes(group = sample), size = 0.5, alpha = 0.05,position = position_jitterdodge(0.3))+
  #ggplot2::geom_jitter(ggplot2::aes(group = sample),width = 0.25,size = 0.5, alpha = 0.05) +
  #geom_jitter(aes(colour = group), width = 0.3, size=2, alpha=0.5) geom_point
  ggplot2::geom_smooth(method = "loess", size = 0.5, se = T) +
  ggplot2::scale_x_continuous(breaks = c(1,2),labels=c("6","13")) +
  ggplot2::ylab(expression(paste("scaled expression"))) +
  ggplot2::xlab("time (dpf)") +
  ggplot2::facet_wrap(~cluster,scales = "free_y",ncol = 3)
tc_clust_plot+theme_classic()+ scale_y_continuous(limits =c(0,2.2))+ 
  scale_color_manual(values = c("violet", "royalblue"))+ 
  theme(legend.title = element_blank())+ 
  theme(legend.position = "none")
dev.off()

#cluster46
tiff(file = "./figures/neurogenesis_cl_46.tiff",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 10, res=1200, units = "cm", compression = "lzw")
 # The height of the plot in inches

#tc_data_long.wt.pos$cluster <- factor(tc_data_long.wt.pos$cluster,      # Reordering group factor levels
#                                      levels = c("cluster 2","cluster 7","cluster 4","cluster 8","cluster 1","cluster 3","cluster 5","cluster 6","cluster 9"))

cluster46= filter(tc_data_long.wt.pos,cluster %in% c("cluster 4","cluster 6"))
tc_clust_plot <- ggplot2::ggplot(cluster46, ggplot2::aes(timepoint, cpm, colour = genotype)) +
  ggplot2::geom_hline(color = "gray", yintercept = 0, linetype = 2, size = 0.5) +
  #ggplot2::geom_line(ggplot2::aes(group = trace), size = 1, alpha = 0.05) +
  ggplot2::geom_boxplot(ggplot2::aes(group = sample), size = 0.5,alpha = 0.1,outlier.size = 0.5) +
  ggplot2::geom_point(ggplot2::aes(group = sample), size = 0.5, alpha = 0.05,position = position_jitterdodge(0.3))+
  #ggplot2::geom_jitter(ggplot2::aes(group = sample),width = 0.25,size = 0.5, alpha = 0.05) +
  #geom_jitter(aes(colour = group), width = 0.3, size=2, alpha=0.5) geom_point
  ggplot2::geom_smooth(method = "loess", size = 0.5, se = T) +
  ggplot2::scale_x_continuous(breaks = c(1,2),labels=c("6","13")) +
  ggplot2::ylab(expression(paste("scaled expression"))) +
  ggplot2::xlab("time (dpf)") +
  ggplot2::facet_wrap(~cluster,scales = "free_y",ncol = 1)
tc_clust_plot+theme_classic()+ scale_y_continuous(limits =c(0,2.2))+ 
  scale_color_manual(values = c("violet", "royalblue"))+ 
  theme(legend.title = element_blank())+ 
  theme(legend.position = "none")
dev.off()

