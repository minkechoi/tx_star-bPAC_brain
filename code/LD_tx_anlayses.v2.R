
# read_tables -------------------------------------------------------------

time.name2 = c("m4.s", "m4_LD.s", "s.m4-LD","t.m4-LD")

for(i in 9:10){
  deg = get(paste0("deg.cpm",i))
  tb = deg
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
  tb[,"source"] = time.name2[i-6]
  name = paste0("deg.cpm",i)
  assign(name,tb)
}

#write.csv(deg.cpm9, "./data/deg.cpm9.csv")
#write.csv(deg.cpm10, "./data/deg.cpm10.csv")

# known GC responding genes -----------------------------------------------

#GC pathway,nCPM
library(plotrix)

GC_path = c("nr3c1","nr3c2","igfbp1a","fkbp5","tsc22d3","hsd11b2")

# gc_path
gcp.deg = gene_all%>% filter(zfin_id_symbol %in% GC_path)%>%arrange(factor(zfin_id_symbol, level = GC_path ))%>% 
  filter(source %in% time.name[c(7:8)])
GC_pathway = c(rep("nr3c1",4),rep("nr3c2",4),rep("igfbp1a",4),rep("fkbp5",4),rep("tsc22d3",4),rep("hsd11b2",4))
time_point <- rep(c(rep("d120",2),rep("LD",2)),6)
geno = rep(c("wt","bPAC+"),12)

mean_CPM = c()
for (i in 0:11) {
  mean_CPM[(2*i)+1]=gcp.deg[i+1,21]
  mean_CPM[(2*i)+2]=gcp.deg[i+1,20]
  
}

se=c()
for (i in 0:11) {
  se[(2*i)+1]=std.error(t(gcp.deg[i+1,c(15:19)]))
  se[(2*i)+2]=std.error(t(gcp.deg[i+1,c(10:14)]))
}

data <- data.frame(GC_pathway,time_point,geno,mean_CPM,se)
data[,"id"]=paste0(time_point,"_",geno)

data$GC_pathway=factor(data$GC_pathway, level = GC_path)
data$time_point=factor(data$time_point, level = c("d120","LD"))
data$geno=factor(data$geno, level = c("wt","bPAC+"))
data$id=factor(data$id, level = c("d120_wt","d120_bPAC+", "LD_wt","LD_bPAC+" ))

##
library(viridis)
require(quantmod)
require(ggplot2)
library(ggsignif)
desturate = rep(c("grey","#FFA0A0"),2)
# Grouped

a= ggplot(data, aes(fill=id, y=mean_CPM, x=time_point )) + 
  #geom_hline(yintercept=c(1,-1),col="gray",linetype="dashed", size=0.5)+
  #geom_hline(yintercept=c(0),col="black",linetype="solid", size=0.5)+
  geom_bar(position="dodge", stat="identity",alpha=0.8, col="black")+
  geom_errorbar(aes(ymin=mean_CPM, ymax=mean_CPM+se), width=0.2,
                position=position_dodge(.9))+
  geom_signif(comparisons = list(c("d120", "LD")), annotations=c("*"), y_position = 18, tip_length = 0.03,textsize = 8) +
  #geom_rect(aes(xmin = 3-.5, xmax = 4+.5,
  #              ymin = -Inf, ymax = Inf,fill = time_point), alpha = .2)+
  theme_classic() +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 15,face="bold", vjust = 0.2),
        legend.position = "none",
        axis.text.x = element_text(size = 13,face="italic"),
        strip.text.x = element_text(size = 15,face = "italic"))+
  scale_fill_manual(values = desturate)+
  scale_y_continuous(expand=expansion(mult = c(0, 0.4)))+
  labs(title = "",x="")+
  facet_wrap(~GC_pathway,scales = "free_y",ncol = 3)


#Error bar correction
gcp.deg.ld = rbind(deg.cpm7,deg.cpm8,deg.cpm9,deg.cpm10)%>% filter(zfin_id_symbol %in% GC_path)%>%arrange(factor(zfin_id_symbol, level = GC_path ))

pval= c(gcp.deg.ld$t.test)
pval_name= c(paste0(gcp.deg.ld$zfin_id_symbol,"_",gcp.deg.ld$source) )
point=gcp.deg.ld$source
pval_as = data.frame(pval_name,pval,point)

pval_as$pval = ifelse(pval_as$pval <= 0.001, "***",
                      ifelse(pval_as$pval <= 0.01, "**",
                             ifelse(pval_as$pval < 0.05, "*", NA)))
pval_as=na.omit(pval_as)

# revise * bars
myplot= ggplot_build(a)
tem.t= myplot$data[[3]]
tem.t=rbind(tem.t,tem.t,tem.t,tem.t[1:3,])
myplot$data[[3]]=tem.t

myplot$data[[3]]$annotation = c(pval_as$pval,pval_as$pval,pval_as$pval)
myplot$data[[3]]$group = c(pval_as$pval_name,pval_as$pval_name,pval_as$pval_name)
myplot$data[[3]]$PANEL = c(rep(c(1,1, #nr3c1
                                 2,2,2, #nr3c2
                                 3,3,3,3, #igfbp1a
                                 4,4,4, # fkbp5
                                 5,5,5,5, # tsc22d3
                                 6,6,6), #hsd11b2
                               3)) # 3 components for each bar

x= pval_as$point
x= gsub(pattern ="t.m4-LD",x = x,replacement = 0.75)
x= gsub(pattern ="s.m4-LD",x = x,replacement = 1.25)
x= gsub(pattern = "d120",x = x,replacement = 0.75)
x= gsub(pattern ="LD",x = x,replacement = 1.75)
x= as.numeric(x)

xend= pval_as$point
xend= gsub(pattern ="t.m4-LD",x = xend,replacement = 0.75+1)
xend= gsub(pattern ="s.m4-LD",x = xend,replacement = 1.25+1)
xend= gsub(pattern = "d120",x = xend,replacement = 0.75+0.5)
xend= gsub(pattern ="LD",x = xend,replacement = 1.75+0.5)
xend= as.numeric(xend)

myplot$data[[3]]$x=c(x,x,xend)
myplot$data[[3]]$xend=c(xend,x,xend)



y=c(150,170, #nr3c1
    57,57,65, #nr3c2
    8,15,19,12, #igfbp1a
    380,600,520, #fkbp5
    230,360,420,300, #tsc22d3
    47.5,52,44) # hsd11b2
myplot$data[[3]]$y = c(y,
                       y[1]-5,y[2]-5, #nr3c1
                       y[3]-2,y[4]-2,y[5]-2, #nr3c2
                       y[6]-0.8,y[7]-0.8,y[8]-0.8,y[9]-0.8, #igfbp1a
                       y[10]-20,y[11]-20,y[12]-20, #fkbp5
                       y[13]-10,y[14]-10,y[15]-10,y[16]-10, #tsc22d3
                       y[17]-(5/3),y[18]-(5/3),y[19]-(5/3),  # hsd11b2
                       y[1]-5,y[2]-5, #nr3c1
                       y[3]-2,y[4]-2,y[5]-2, #nr3c2
                       y[6]-0.8,y[7]-0.8,y[8]-0.8,y[9]-0.8, #igfbp1a
                       y[10]-20,y[11]-20,y[12]-20, #fkbp5
                       y[13]-10,y[14]-10,y[15]-10,y[16]-10, #tsc22d3
                       y[17]-(5/3),y[18]-(5/3),y[19]-(5/3)  # hsd11b2
                       ) 

myplot$data[[3]]$yend = c(y,y,y)

myplot$data[[3]]$textsize =c(rep(8,nrow(myplot$data[[3]])))
myplot$data[[3]]$vjust =c(rep(0.5,nrow(myplot$data[[3]])))
#myplot$data[[3]]$colour[c(5,5+19,5+19+19,19,19*2,19*3)]="red"
a=ggplot_gtable(myplot)
plot(a)


tiff("./figures/Fig4C_GC_genes.24x15.tiff", width = 20,height =15, units = "cm", res=1200, pointsize = 4, compression = "lzw" )
plot(a)
dev.off() 



# QC, filtering ---------------------------------------------------------------


for(i in 9:10){
  deg=get(paste0("deg.cpm",i))
  
  deg = deg %>% filter(mean.cpm.pos > 5|mean.cpm.wt >5)%>% filter(t.test < 0.05)%>%
    filter(FDR < 0.05)%>% arrange(desc(logCPM))
  d = duplicated(deg$zfin_id_symbol)
  deg=deg[!d,]
  name = paste0("ft.deg.cpm",i)
  assign(name,deg)
}


#LD_DEG calling -------------------------------------------------------------

#LD.DEGs

for (i in 9:10) {
  deg=get(paste0("ft.deg.cpm",i))
  deg.up = deg %>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)%>% filter(logFC > log2(1.5)) %>% 
    filter(FDR < 0.01)%>% arrange(desc(logFC))
  deg.dn = deg %>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)%>% filter(logFC < -log2(1.5))  %>% 
    filter(FDR < 0.01)%>% arrange(logFC)
  name1 = paste0("ft.deg.cpm.up",i)
  assign(name1,deg.up)
  name2 = paste0("ft.deg.cpm.dn",i)
  assign(name2,deg.dn)
  
}

### isolation of cpm values of wt and star:bPAC+

for(i in c(5:8)){
  deg=get(paste0("deg.cpm",i))
  deg= filter(deg, id %in% b.ft.id)
  deg = deg[,c(8,10:19)]
  name = paste0("wt.pos.cpm",i)
  assign(name,deg)
}

wt.pos.cpm = wt.pos.cpm5

for(i in c(6:8)){
  deg=get(paste0("wt.pos.cpm",i))
  wt.pos.cpm = cbind(wt.pos.cpm,deg[,c(2:11)])
  
}

colnames(wt.pos.cpm)= c("zfin_id_id",
                        paste0(time.name[6],"_","pos",1:5),
                        paste0(time.name[6],"_","wt",1:5),
                        paste0(time.name[7],"_","pos",1:5),
                        paste0(time.name[7],"_","wt",1:5),
                        paste0(time.name[9],"_","pos",1:5),
                        paste0(time.name[9],"_","wt",1:5),
                        paste0(time.name[10],"_","pos",1:5),
                        paste0(time.name[10],"_","wt",1:5)
)

wt.pos.cpm = na.omit(wt.pos.cpm)
#write.csv(wt.pos.cpm, "./data/LD_wt.pos.cpm.csv")

#selection of DEG established in early life
sig1 = filter(deg.cpm5, FDR < 0.01) %>% filter(mean.cpm.pos > 5|mean.cpm.wt >5)%>% filter(t.test < 0.05)%>%
  filter(abs(logFC) >= log2(2))%>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)
sig2 = filter(deg.cpm6, FDR < 0.01) %>% filter(mean.cpm.pos > 5|mean.cpm.wt >5)%>% filter(t.test < 0.05)%>%
  filter(abs(logFC) >= log2(2))%>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)
sig3 = filter(deg.cpm7, FDR < 0.05)

#DEG after LD between star:bPAC and wt
sig4 = deg.cpm8 %>% filter(FDR < 0.01) %>% filter(mean.cpm.pos > 5|mean.cpm.wt >5)%>% filter(t.test < 0.05)%>%
  filter(abs(logFC) >= log2(2))%>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)

#DEG after LD in star:bPAC
sig5 = filter(deg.cpm9, FDR < 0.01)%>% filter(mean.cpm.pos > 5|mean.cpm.wt >5)%>% filter(t.test < 0.05)%>%
  filter(abs(logFC) >=  log2(1.5))%>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)

#non-de after LD in star:bPAC
sig5.non = filter(deg.cpm9, FDR > 0.01) %>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)%>% filter(mean.cpm.pos > 5|mean.cpm.wt >5)%>% filter(t.test < 0.05)

#DEG after LD in TU
sig6 = filter(deg.cpm10, FDR < 0.01)%>% filter(mean.cpm.pos > 5|mean.cpm.wt >5)%>% filter(t.test < 0.05)%>%
  filter(abs(logFC) >=  log2(1.5))%>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)
#Common DEGs following LD
#overlapping tu vs.tuLD between bPAC vs. bPAC+LD
sig56_genes=intersect(sig5$zfin_id_symbol,sig6$zfin_id_symbol)#44

#check direction of FC
sig5$zfin_id_symbol[c(which(sig5[c(which(sig5$zfin_id_symbol%in%sig56_genes)),"logFC"]*sig6[c(which(sig6$zfin_id_symbol%in%sig56_genes)),"logFC"] <0))]
#if it's 0. proceed to next

sig56_5 = sig5[c(which(sig5$zfin_id_symbol%in%sig56_genes)),]
sig56_6 = sig6[c(which(sig6$zfin_id_symbol%in%sig56_genes)),]
sig56 = rbind(sig56_5,sig56_6)

#check overlapping with LD
deg.cpm8[c(which(deg.cpm8$zfin_id_symbol%in%sig56_genes)),] %>% 
  filter(abs(logFC)>log2(2)) %>% filter(FDR <0.01)
#if it's 0. proceed to next

  #.. overlapping 44 genes show same direction of regulation, and are not part of post-LD-degs 
  #.. This result indicates that those genes are involved in normal response to LD.

deg12= unique(c(deg.list$`DEG 2`,deg.list$`DEG 1`))
deg12.primed = deg12[-c(which(deg12 %in% c(intersect(deg.list$`DEG 2`,deg.list$`DEG 1`))))]


#check direction of FC
sig46_genes=intersect(sig4$zfin_id_symbol,sig6$zfin_id_symbol) # "csrnp1b" "lipg"
sig6$zfin_id_symbol[c(which(sig6[c(which(sig6$zfin_id_symbol%in%sig46_genes)),"logFC"]*sig6[c(which(sig4$zfin_id_symbol%in%sig46_genes)),"logFC"] <0))]
  #.. DEG1 and DEG3 have 2 overlapping gene but oppositely regulated. Thus, they are part of post-LD_DEGs

ab.LD.res = unique(c(sig5$id,sig4$id,sig6$id))
ab.LD.res=ab.LD.res[- c(which(ab.LD.res %in% unique(sig56$id)))]
ab.LD.res.deg.gene = intersect(deg.cpm5[c(which(deg.cpm5$id %in% ab.LD.res )),8],prot_id_map$SYMBOL)

LD.deg.list = rbind(sig4,sig5,sig6)
LD.deg.list = LD.deg.list[-c(which(LD.deg.list$zfin_id_symbol %in% unique(sig56$zfin_id_symbol))),]
LD.deg.list$source= gsub("s.m4-LD",x = LD.deg.list$source ,"DEG2")
LD.deg.list$source= gsub("t.m4-LD",x = LD.deg.list$source ,"DEG1")
LD.deg.list$source= gsub("LD",x = LD.deg.list$source ,"DEG3")


#s.table3_list of LD-DEGs
write.csv(LD.deg.list, "./outputs/s.table3_list of LD-DEGs.csv")
####
#####

#generating tables for heatmap

for (i in c(7,8)) {
  deg=get(paste0("deg.cpm",i))
  deg = deg %>% filter(zfin_id_symbol %in% ab.LD.res.deg.gene)%>% arrange(desc(zfin_id_symbol))
  md.fc = deg[,c("mean.cpm.pos", "mean.cpm.wt")]
  rownames(md.fc)=deg$zfin_id_symbol
  name1 = paste0("ad.heat.deg.cpm",i)
  assign(name1,md.fc)
}
ht.mt.ad = data.frame(ad.heat.deg.cpm7,
                      ad.heat.deg.cpm8)
colnames(ht.mt.ad)=c("120_bPAC+","120_wt","LD_bPAC+","LD_wt")

st_sig_tx =as.matrix(t(scale(t(ht.mt.ad))))
rownames(st_sig_tx)=rownames(ht.mt.ad)

distfunc <- function(x) dist(x, method="euclidean")
hclustfunc <- function(x) hclust(x, method="ward.D2")

hc = hclust(dist(t(st_sig_tx), method="euclidean"), method="ward.D2")
dd = as.dendrogram(hc)

#fviz_nbclust(st_sig_tx, FUN = hcut, method = "silhouette")

rhc = hclust(dist(st_sig_tx, method="euclidean"), method="ward.D2")
gr.row <- cutree(rhc, 2) 

col1 <- brewer.pal(3, "Set3")
desturate =viridis(4)

cpm.table = data.frame(rhc$labels,rhc$order,gr.row)
#write.csv(cpm.table,"./outputs/cpm.table_heat_4clussters.csv")

png("./figures/Fig5B_ad.tree_heatmap.png", width = 300, height = 200)
plot(dd)
dev.off()

tiff("./figures/Fig4B.ad_heatmap.6x12.tiff", width = 6, height = 12, res = 300, units = "cm", compression = "lzw")
hvld = heatmap.2((st_sig_tx), col=mycol, scale='none',
                 dendrogram = "both",breaks=breaks,hclustfun=hclustfunc, distfun=distfunc,
                 #main=c(paste0("Z-score (",nrow(st_sig_tx),")")), Colv = reorder(as.dendrogram(dd),4:1),
                 density='none', trace = 'none',
                 na.color = F,
                 offsetRow = 0,
                 offsetCol = c(0,0),
                 labRow = "",
                 colsep = c(2),
                 cexCol = 1,
                 #labCol = c("13_pos","6_pos","6_wt","13_wt","LD_wt","120_wt","120_pos","LD_pos"),
                 #cellnote=ht.FDR,
                 #notecex=1,
                 #notecol="black",
                 ColSideColors = c(rep(c(desturate[1],desturate[4]),2)),
                 RowSideColors=col1[gr.row],
                 margins = c(5, 5),
                 key.title = "",
                 key.xlab = "",
                 key.ylab = "",
                 keysize = 0.5,
                 key = F)

dev.off()

nrow(st_sig_tx)


#### venn diagram
deg.list = list(
  "DEG 1" = sig6$zfin_id_symbol,
  "DEG 2" = sig5$zfin_id_symbol,
  "DEG 3" = sig4$zfin_id_symbol)


desturate=viridis(3)

library(eulerr)


fit2= euler(deg.list,shape = "ellipse")
fit2$fitted.values[5]=2
fit2$ellipses$k[1]=-1.1
fit2$ellipses$a[1]=1.7
fit2$ellipses$b[1]=6
fit2$ellipses$phi[1]=3.1
deg.venn =  plot(fit2, 
      quantities = F,
      fill = viridis(3), alpha =0.3,
      lty = 1, 
      labels = ""
)

#




tiff("./figures/Fig5E_venn.LD.deg.12x12.tiff",  height = 8,width = 8, res=1200, units = "cm", compression = "lzw")
plot(deg.venn)
dev.off()

'                 original fitted residuals regionError
DEG 1                   18     18         0       0.000
DEG 2                  463    463         0       0.000
DEG 3                  922    922         0       0.000
DEG 1&DEG 2             44     44         0       0.000
DEG 1&DEG 3              2      2         2       0.001
DEG 2&DEG 3           1110   1110         0       0.000
DEG 1&DEG 2&DEG 3        0      0         0       0.000

diagError: 0.001 
stress:    0 
'

#####
#GC-primed genes
deg12= unique(c(deg.list$`DEG 2`,deg.list$`DEG 1`))
deg12.primed = deg12[-c(which(deg12 %in% c(intersect(deg.list$`DEG 2`,deg.list$`DEG 1`))))]


### not working well. thus, used below
"Luana Micallef and Peter Rodgers (2014). eulerAPE: Drawing Area-proportional 3-Venn Diagrams Using Ellipses.
PLoS ONE 9(7): e101717. doi:10.1371/journal.pone.0101717. http://www.eulerdiagrams.org/eulerAPE
"

# chea3_TF enrichment test ------------------------------------------------

###chea3_results
#Chea3 API: https://maayanlab.cloud/chea3/
# run ChIP-X Enrichment analysis version 3

library(httr)
library(jsonlite)
#up-regulated
LD.deg_all.up=LD.deg.list%>%filter(logFC > 0)%>%filter(zfin_id_symbol %in% deg12.primed)
zh.LD_DEG.up= unique(filter(zh.all, zfin_id_symbol %in% LD.deg_all.up$zfin_id_symbol)[,2])

url = "https://maayanlab.cloud/chea3/api/enrich/"
encode = "json"
payload = list(query_name = "LD_DEG_up", gene_set =zh.LD_DEG.up )

#POST to ChEA3 server
response = POST(url = url, body = payload, encode = encode)
json = content(response, "text")  

chea3_LD = fromJSON(json)
overlapped.number= str_count(chea3_LD$`Integrated--meanRank`$Overlapping_Genes,",")+1
chea3_LD.int.mean.rank = cbind(chea3_LD$`Integrated--meanRank`,overlapped.number)
write.csv(chea3_LD.int.mean.rank,"./outputs/s.table5_chea3_LD_DEG.up_results_Integrated--meanRank.csv")  


chea3 = data.frame(chea3_LD$`Integrated--meanRank`)
chea3.sb50= chea3 %>% filter(as.numeric(Rank) <= 50) 
chea3.up.deg = chea3 %>% filter(as.numeric(Rank) <= 50) %>%
  filter(TF %in% c(zh.all$human.Symbol[c(which(zh.all$zfin_id_symbol %in% deg_all.up$zfin_id_symbol))]))

library(forcats) 

top.tf.LD=chea3.up.deg$TF
top.tf.LD.rank= paste0("rank:",chea3.up.deg$Rank)
top.tf.LD.score= log2(1632/as.numeric(chea3.up.deg$Score))

#.. is part of LD_DEG? 
topTF.zh_LD=filter(LD.deg.list, zfin_id_symbol %in% unique(c(filter(zh.all,human.Symbol %in% top.tf.LD)$zfin_id_symbol)))$zfin_id_symbol
topTF.hz_LD=unique(c(filter(zh.all,zfin_id_symbol %in% topTF.zh_LD)$human.Symbol))

#.. "CSRNP1" "EGR1"   "FOSL2"  "JUNB"   "EN2"    "FEZF2"  "FOS"   ,for up-reg LD-dges.  
#.. "FEZF2"  "CSRNP1" "FOS"    "BATF2"  are part of LD_DEGs

#top10 TF .but only 9 are detected
tf.bp.data = data.frame("TF"=top.tf.LD[1:7],"rank"=top.tf.LD.rank[1:7],"score"=top.tf.LD.score[1:7])

col= c()
for (i in 1:7) {
  if (tf.bp.data$TF[i] %in% topTF.hz_LD) {
    col= c(col,"forestgreen")
  }else {
    col= c(col,"grey")
  }
}
tf.bp.data[,"col"]=col
tf.bp.data=mutate(tf.bp.data,TF = fct_reorder(TF, score))

##barplot

LD.tf.p= ggplot(tf.bp.data, aes(x=TF, y=score))+
  geom_bar(stat = "identity",fill=col, alpha=0.5, width = 0.9 )+
  geom_text(aes(label = rank), stat = "identity", hjust = 1.5, colour = "black")+
  coord_flip()+
  xlab("")+
  ylab("TF enrichment score")+
  theme_classic() +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20,face=c("italic")),
        axis.text.x = element_text(size = 15),
        #axis.text.y = element_blank(),
        axis.title.x = element_text(size = 15),
        strip.text.x = element_text(face = "italic"),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  #labs(title = "",x="",y="cortisol (pg/larva)")+
  #scale_fill_manual(values = desturate)+
  scale_y_continuous(expand=expansion(mult = c(0, 0.3)))


tiff("./figures/s.Fig6a_LD-up.tf.tiff", width = 12, height = 12,units = "cm", res = 1200, compression = "lzw")
plot(LD.tf.p)
dev.off()

#down-regulated
LD.deg_all.dn=LD.deg.list%>%filter(logFC < 0)%>%filter(zfin_id_symbol %in% deg12.primed)

zh.LD_DEG.dn= unique(filter(zh.all, zfin_id_symbol %in% LD.deg_all.dn$zfin_id_symbol)[,2])


url = "https://maayanlab.cloud/chea3/api/enrich/"
encode = "json"
payload = list(query_name = "down_LD_DEG", gene_set =zh.LD_DEG.dn )

#POST to ChEA3 server
response = POST(url = url, body = payload, encode = encode)
json = content(response, "text")  

chea3_LD = fromJSON(json)
overlapped.number= str_count(chea3_LD$`Integrated--meanRank`$Overlapping_Genes,",")+1
chea3_LD.int.mean.rank = cbind(chea3_LD$`Integrated--meanRank`,overlapped.number)
write.csv(chea3_LD.int.mean.rank,"./outputs/s.table5_chea3_LD_DEG.dn_results_Integrated--meanRank.csv")  

chea3 = data.frame(chea3_LD$`Integrated--meanRank`)
chea3.sb50= chea3 %>% filter(as.numeric(Rank) <= 50) 
chea3.dn.deg = chea3 %>% filter(as.numeric(Rank) <= 50) %>%
  filter(TF %in% c(zh.all$human.Symbol[c(which(zh.all$zfin_id_symbol %in% deg_all.dn$zfin_id_symbol))]))

library(forcats) 

top.tf.LD=chea3.dn.deg$TF
top.tf.LD.rank= paste0("rank:",chea3.dn.deg$Rank)
top.tf.LD.score= log2(1632/as.numeric(chea3.dn.deg$Score))

#.. is part of LD_DEG? 
topTF.zh_LD=filter(LD.deg.list, zfin_id_symbol %in% unique(c(filter(zh.all,human.Symbol %in% top.tf.LD)$zfin_id_symbol)))$zfin_id_symbol
topTF.hz_LD=unique(c(filter(zh.all,zfin_id_symbol %in% topTF.zh_LD)$human.Symbol))

#top10 TF
tf.bp.data = data.frame("TF"=top.tf.LD[1:10],"rank"=top.tf.LD.rank[1:10],"score"=top.tf.LD.score[1:10])

col= c()
for (i in 1:10) {
  if (tf.bp.data$TF[i] %in% topTF.hz_LD) {
    col= c(col,"forestgreen")
  }else {
    col= c(col,"grey")
  }
}
tf.bp.data[,"col"]=col
tf.bp.data=mutate(tf.bp.data,TF = fct_reorder(TF, score))

##barplot

LD.tf.p= ggplot(tf.bp.data, aes(x=TF, y=score))+
  geom_bar(stat = "identity",fill=col, alpha=0.5, width = 0.9 )+
  geom_text(aes(label = rank), stat = "identity", hjust = 1.5, colour = "black")+
  coord_flip()+
  xlab("")+
  ylab("TF enrichment score")+
  theme_classic() +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20,face=c("italic")),
        axis.text.x = element_text(size = 15),
        #axis.text.y = element_blank(),
        axis.title.x = element_text(size = 15),
        strip.text.x = element_text(face = "italic"),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  #labs(title = "",x="",y="cortisol (pg/larva)")+
  #scale_fill_manual(values = desturate)+
  scale_y_continuous(expand=expansion(mult = c(0, 0.3)))


tiff("./figures/s.Fig6b_LD-down.tf.tiff", width = 12, height = 12,units = "cm", res = 1200, compression = "lzw")
plot(LD.tf.p)
dev.off()



# GO enrichment analysis of LD-DEGs ---------------------------------------


#list of background genes
for (i in c(8:10)){
  deg = get(paste0("deg.cpm",i))
  deg = deg %>% filter(mean.cpm.pos > 5|mean.cpm.wt >5)
  deg = prot_id_map[deg$zfin_id_symbol, "ENTREZID"]
  deg = na.omit(deg)
  name = paste0("bk.gene.LD",i)
  assign(name, deg)
  print(i)
}

bk.all.LD = unique(c(bk.gene.LD8,bk.gene.LD9,bk.gene.LD10))

#each time point
time.name2=c("DEG","DEG1","DEG2","DEG3")

#up
tc_point_up <- list()
## Reformat the timecourse cluster classifications as a list

for (i in c(2:4)){
  cluster_psites <- which(LD.deg_all.up$source == time.name2[i])
  cluster_prots <- unique(LD.deg_all.up$zfin_id_symbol[cluster_psites])
  cluster_prots_entrez <- prot_id_map[cluster_prots, "ENTREZID"]
  cluster_prots_entrez=na.omit(cluster_prots_entrez)
  tc_point_up[[time.name2[i]]] <- cluster_prots_entrez
}

#down
tc_point_dn <- list()
## Reformat the timecourse cluster classifications as a list


for (i in c(2:4)){
  cluster_psites <- which(LD.deg_all.dn$source == time.name2[i])
  cluster_prots <- unique(LD.deg_all.dn$zfin_id_symbol[cluster_psites])
  cluster_prots_entrez <- prot_id_map[cluster_prots, "ENTREZID"]
  cluster_prots_entrez=na.omit(cluster_prots_entrez)
  tc_point_dn[[time.name2[i]]] <- cluster_prots_entrez
}


LD.tc_enrich.up.tb = data.frame()
LD.tc_enrich.dn.tb = data.frame()

gos=c("BP","MF","CC")


for (k in 1:3) {
  ##up
  tryCatch({
    tc_enrich.up <- clusterProfiler::compareCluster(tc_point_up, fun = "enrichGO",  OrgDb='org.Dr.eg.db', ont = gos[k],
                                                  bk.all.LD,pvalueCutoff = 0.05,
                                                  pAdjustMethod = "BH",qvalueCutoff = 0.2,
                                                  minGSSize = 10,
                                                  maxGSSize = 600, 
                                                  readable = TRUE)
           
  tc_enrich.up.tb = data.frame(tc_enrich.up[1:length(tc_enrich.up[]$ID)])
  tc_enrich.up.tb[,"source"]=gos[k]
  
  a= str_split(tc_enrich.up.tb[,4],pattern = "/")
  if (length(a)>1) {
    b= data.frame(Reduce(rbind,a))
    } else {b= data.frame(a[[1]][1],a[[1]][2])  
  }
  
  tc_enrich.up.tb[,"n.deg"] = as.numeric(b[,2])
  
  a= str_split(tc_enrich.up.tb[,5],pattern = "/")
  if (length(a)>1) {
    b= data.frame(Reduce(rbind,a))
  } else {b= data.frame(a[[1]][1],a[[1]][2])  
  }
  tc_enrich.up.tb[,"gene.in.catg"] = as.numeric(b[,1])
  
  tc_enrich.up.tb[,"ratio.catg"] = tc_enrich.up.tb[,10]/tc_enrich.up.tb[,"gene.in.catg"]
  
  #filtering
  tc_enrich.up.tb= tc_enrich.up.tb %>% filter(Count >= 10) %>% filter(ratio.catg >= 0.05)
  
  LD.tc_enrich.up.tb =rbind(LD.tc_enrich.up.tb,tc_enrich.up.tb)
  
  
  ##down
  tc_enrich.dn <- clusterProfiler::compareCluster(tc_point_dn, fun = "enrichGO",  OrgDb='org.Dr.eg.db', ont = gos[k],
                                                  bk.all.LD,pvalueCutoff = 0.05,
                                                  pAdjustMethod = "BH",qvalueCutoff = 0.2,
                                                  minGSSize = 10,
                                                  maxGSSize = 1000, 
                                                  readable = TRUE)
           
  tc_enrich.dn.tb = data.frame(tc_enrich.dn[1:length(tc_enrich.dn[]$ID)])
  tc_enrich.dn.tb[,"source"]=gos[k]
  
  a= str_split(tc_enrich.dn.tb[,4],pattern = "/")
  if (length(a)>1) {
    b= data.frame(Reduce(rbind,a))
  } else {b= data.frame(a[[1]][1],a[[1]][2])  
  }
  tc_enrich.dn.tb[,"n.deg"] = as.numeric(b[,2])
  
  a= str_split(tc_enrich.dn.tb[,5],pattern = "/")
  if (length(a)>1) {
    b= data.frame(Reduce(rbind,a))
  } else {b= data.frame(a[[1]][1],a[[1]][2])  
  }
  tc_enrich.dn.tb[,"gene.in.catg"] = as.numeric(b[,1])
  
  tc_enrich.dn.tb[,"ratio.catg"] = c(tc_enrich.dn.tb[,10]/tc_enrich.dn.tb[,"gene.in.catg"])
  
  #filtering
  tc_enrich.dn.tb= tc_enrich.dn.tb %>% filter(Count >= 10) %>% filter(ratio.catg >= 0.05)
  
  LD.tc_enrich.dn.tb =rbind(LD.tc_enrich.dn.tb, tc_enrich.dn.tb)}, error=function(e){})
}


#write_combined_table

write.csv(LD.tc_enrich.up.tb,"./outputs/s.table4a_LD.tc_enrich.up.tb.csv" )
write.csv(LD.tc_enrich.dn.tb,"./outputs/s.table4b_LD.tc_enrich.dn.tb.csv" )



#all_times LD-DEGs
#up

tc_point_up <- list()
## Reformat the timecourse cluster classifications as a list

  cluster_prots <- unique(LD.deg_all.up$zfin_id_symbol)
  cluster_prots_entrez <- prot_id_map[cluster_prots, "ENTREZID"]
  cluster_prots_entrez=na.omit(cluster_prots_entrez)
  tc_point_up[["up"]] <- cluster_prots_entrez


#down
tc_point_dn <- list()
## Reformat the timecourse cluster classifications as a list

  cluster_prots <- unique(LD.deg_all.dn$zfin_id_symbol)
  cluster_prots_entrez <- prot_id_map[cluster_prots, "ENTREZID"]
  cluster_prots_entrez=na.omit(cluster_prots_entrez)
  tc_point_dn[["down"]] <- cluster_prots_entrez



all.LD.tc_enrich.up.tb = data.frame()
all.LD.tc_enrich.dn.tb = data.frame()
gos=c("BP","MF","CC")

#warning, duplication in the result if there is no result. 

for (k in 1:3) {
  tryCatch({
  ##up
  tc_enrich.up <- clusterProfiler::compareCluster(tc_point_up, fun = "enrichGO",  OrgDb='org.Dr.eg.db', ont = gos[k],
                                                  bk.all.LD,pvalueCutoff = 0.05,
                                                  pAdjustMethod = "BH",qvalueCutoff = 0.2,
                                                  minGSSize = 10,
                                                  maxGSSize = 600, 
                                                  readable = TRUE)
  
  tc_enrich.up.tb = data.frame(tc_enrich.up[1:length(tc_enrich.up[]$ID)])
  tc_enrich.up.tb[,"source"]<-gos[k]
  
  a= str_split(tc_enrich.up.tb[,4],pattern = "/")
  if (length(a)>1) {
    b= data.frame(Reduce(rbind,a))
  } else {b= data.frame(a[[1]][1],a[[1]][2])  
  }  
  tc_enrich.up.tb[,"n.deg"] = as.numeric(b[,2])
  
  a= str_split(tc_enrich.up.tb[,5],pattern = "/")
  if (length(a)>1) {
    b= data.frame(Reduce(rbind,a))
  } else {b= data.frame(a[[1]][1],a[[1]][2])  
  }  
  tc_enrich.up.tb[,"gene.in.catg"] = as.numeric(b[,1])
  
  tc_enrich.up.tb[,"ratio.catg"] = tc_enrich.up.tb[,10]/tc_enrich.up.tb[,"gene.in.catg"]
  
  #filtering
  tc_enrich.up.tb= tc_enrich.up.tb %>% filter(Count >= 10) %>% filter(ratio.catg >= 0.05)
  
  all.LD.tc_enrich.up.tb =rbind(all.LD.tc_enrich.up.tb,tc_enrich.up.tb)
  
  
  ##down
  tc_enrich.dn <- clusterProfiler::compareCluster(tc_point_dn, fun = "enrichGO",  OrgDb='org.Dr.eg.db', ont = gos[k],
                                                  bk.all.LD,pvalueCutoff = 0.05,
                                                  pAdjustMethod = "BH",qvalueCutoff = 0.2,
                                                  minGSSize = 10,
                                                  maxGSSize = 1000, 
                                                  readable = TRUE)
  
  tc_enrich.dn.tb = data.frame(tc_enrich.dn[1:length(tc_enrich.dn[]$ID)])
  tc_enrich.dn.tb[,"source"]=gos[k]
  
  a= str_split(tc_enrich.dn.tb[,4],pattern = "/")
  if (length(a)>1) {
    b= data.frame(Reduce(rbind,a))
  } else {b= data.frame(a[[1]][1],a[[1]][2])  
  }
  tc_enrich.dn.tb[,"n.deg"] = as.numeric(b[,2])
  
  a= str_split(tc_enrich.dn.tb[,5],pattern = "/")
  if (length(a)>1) {
    b= data.frame(Reduce(rbind,a))
  } else {b= data.frame(a[[1]][1],a[[1]][2])  
  }
  tc_enrich.dn.tb[,"gene.in.catg"] = as.numeric(b[,1])
  
  tc_enrich.dn.tb[,"ratio.catg"] = c(tc_enrich.dn.tb[,10]/tc_enrich.dn.tb[,"gene.in.catg"])
  
  #filtering
  tc_enrich.dn.tb= tc_enrich.dn.tb %>% filter(Count >= 10) %>% filter(ratio.catg >= 0.05)
  
  all.LD.tc_enrich.dn.tb =rbind(all.LD.tc_enrich.dn.tb, tc_enrich.dn.tb) } , error=function(e){})
  
}

#write_combined_table

write.csv(all.LD.tc_enrich.up.tb,"./outputs/s.table4c_all_LD.tc_enrich.up.tb.csv" )
write.csv(all.LD.tc_enrich.dn.tb,"./outputs/s.table4d_all_LD.tc_enrich.dn.tb.csv" )


######sb_LD

for(i in 1:8){
  deg=get(paste0("sb.deg.cpm",i))
  deg = deg %>% filter(zfin_id_symbol %in% sb.613120ld) 
  name = paste0("d6told.sb.deg.cpm",i)
  assign(name,deg)
}

sig.any= unique(c(sb.deg.5$zfin_id_symbol,sb.deg.6$zfin_id_symbol,sb.deg.120$zfin_id_symbol,sb.deg.LD$zfin_id_symbol))


###d6, 13, 120 dpf and LD wt, star:bPAC- and star.bPAC+

ht.mt.allsb = data.frame(d6told.sb.deg.cpm5$logFC,
                         d6told.sb.deg.cpm6$logFC,
                         d6told.sb.deg.cpm7$logFC,
                         d6told.sb.deg.cpm8$logFC,
                         d6told.sb.deg.cpm1$logFC,
                         d6told.sb.deg.cpm2$logFC,
                         d6told.sb.deg.cpm3$logFC,
                         d6told.sb.deg.cpm4$logFC)
row.names(ht.mt.allsb)=d6told.sb.deg.cpm5$zfin_id_symbol

ht.FDR.allsb= data.frame(d6told.sb.deg.cpm5$FDR,
                         d6told.sb.deg.cpm6$FDR,
                         d6told.sb.deg.cpm7$FDR,
                         d6told.sb.deg.cpm8$FDR,
                         d6told.sb.deg.cpm1$FDR,
                         d6told.sb.deg.cpm2$FDR,
                         d6told.sb.deg.cpm3$FDR,
                         d6told.sb.deg.cpm4$FDR)
row.names(ht.FDR.allsb)=d6told.sb.deg.cpm5$zfin_id_symbol

ht.FDR.allsb[ht.FDR.allsb <= 0.05] <-"*"
ht.FDR.allsb[ht.FDR.allsb > 0.05] <-NA
ht.FDR.allsb[abs(ht.mt.allsb) < log2(1.5)] <-NA

st_sig_tx_1= as.matrix(ht.mt.allsb)

distfunc <- function(x) dist(x, method="euclidean")
hclustfunc <- function(x) hclust(x, method="ward.D2")

hc = hclust(dist(t(st_sig_tx_1), method="euclidean"), method="complete")
dd = as.dendrogram(hc)

#fviz_nbclust(st_sig_tx_1, FUN = hcut, method = "silhouette")

rhc = hclust(dist(st_sig_tx_1, method="euclidean"), method="ward.D2")
gr.row <- cutree(rhc, 4) 

col1 <- brewer.pal(4, "Set3")
sb.cluster.table = data.frame(rhc$labels,rhc$order,gr.row)



##### sig.any int. ld-deg
LD_deg= unique(c(LD.deg_all.up$zfin_id_symbol, LD.deg_all.dn$zfin_id_symbol))

sb.int.ld = intersect(LD_deg,sig.any)
sb.rowname = row.names(ht.mt.allsb)
sb.rowname[sb.rowname %in% sb.int.ld] <- "red"
sb.rowname[sb.rowname != "red"] <- "black"
ldrowcols=sb.rowname


tiff("./figures/s.Fig4D_sb.6to120LD.sp.sn.8x30a.tiff", height = 28, width = 8, res = 1200, compression = "lzw", units = "cm")
try(heatmap.2((st_sig_tx_1), col=mycol, scale='none',
               dendrogram = "row",Colv = F,breaks=breaks,hclustfun=hclustfunc, distfun=distfunc,
               #main=c("oxt_regulation, log2FC"),
               density='density', trace = 'none',
               na.color = F,
               offsetRow = 0,
               offsetCol = c(0,0.2),
               cexCol = 1,
               cexRow = 0.9,
               colsep = 4,
               labCol = c("d6+","d13+","d120+", "LD+", "d6-","d13-","d120-", "LD-"),
               #labRow = "",
               colRow = ldrowcols,
               cellnote=ht.FDR.allsb,
               RowSideColors=col1[gr.row],
               notecex=1,
               notecol="black",
               margins = c(5, 5),
               #key.title = "log2FC",
               keysize = 0.5),silent = T)
dev.off()



######Disgenet
# disease associated_enrichment_disignet ----------------------------------
library(disgenet2r)
#disgenet_api_key <- get_disgenet_api_key(
#  email = "user@email.com", 
#  password = "myspwd" )
#Sys.setenv(DISGENET_API_KEY= disgenet_api_key)



#

# for GC-altered
zh_LD_list.deg3= unique(filter(zh.all, zfin_id_symbol %in% c(deg.list$`DEG 3`))[,2])

try(res_enrich <-disease_enrichment( entities =zh_LD_list.deg3, vocabulary = "HGNC",
                                     database = "CURATED" ) # "DEG3"
    ,silent =T)

table.dis.LD <- res_enrich@qresult
table.dis.LD.sig = table.dis.LD %>% arrange(desc(gg)) %>% filter(FDR <0.05) %>% filter(Count >= 20)

a= str_split(table.dis.LD.sig[,4],pattern = "/")
if (length(a)>1) {
  b= data.frame(Reduce(rbind,a))
} else {b= data.frame(a[[1]][1],a[[1]][2])  
}

table.dis.LD.sig[,"n.input"] = as.numeric(b[,2])

a= str_split(table.dis.LD.sig[,5],pattern = "/")
if (length(a)>1) {
  b= data.frame(Reduce(rbind,a))
} else {b= data.frame(a[[1]][1],a[[1]][2])  
}
table.dis.LD.sig[,"gene.in.catg"] = as.numeric(b[,1])

table.dis.LD.sig[,"ratio.catg"] = table.dis.LD.sig[,"Count"]/table.dis.LD.sig[,"gene.in.catg"]

write.csv(table.dis.LD.sig,"./outputs/s.table7_DisGeNET_enrichment_DEG3.csv")

tiff("./figures/Fig6A_top.disgenet.DEG3.tiff", width = 13, height = 12,units = "cm", res = 1200, compression = "lzw")
plot(res_enrich, class = "Enrichment", count =30,  cutoff= 0.05, nchars=70)
dev.off()


# for GC-primed genes
zh_LD_list= unique(filter(zh.all, zfin_id_symbol %in% c(LD.deg_all.dn$zfin_id_symbol,LD.deg_all.up$zfin_id_symbol))[,2])

try(res_enrich <-disease_enrichment( entities =zh_LD_list, vocabulary = "HGNC",
                                 database = "CURATED" ) # "DEG1+2"
,silent =T)
table.dis.LD <- res_enrich@qresult
table.dis.LD.sig = table.dis.LD %>% arrange(desc(gg)) %>% filter(FDR <0.05) %>% filter(Count >= 20)

a= str_split(table.dis.LD.sig[,4],pattern = "/")
if (length(a)>1) {
  b= data.frame(Reduce(rbind,a))
} else {b= data.frame(a[[1]][1],a[[1]][2])  
}

table.dis.LD.sig[,"n.input"] = as.numeric(b[,2])

a= str_split(table.dis.LD.sig[,5],pattern = "/")
if (length(a)>1) {
  b= data.frame(Reduce(rbind,a))
} else {b= data.frame(a[[1]][1],a[[1]][2])  
}
table.dis.LD.sig[,"gene.in.catg"] = as.numeric(b[,1])

table.dis.LD.sig[,"ratio.catg"] = table.dis.LD.sig[,"Count"]/table.dis.LD.sig[,"gene.in.catg"]

write.csv(table.dis.LD.sig,"./outputs/s.table6_DisGeNET_enrichment_DEG12.csv")

tiff("./figures/Fig6A_top.disgenet_DEG12.tiff", width = 13, height = 12,units = "cm", res = 1200, compression = "lzw")
plot(res_enrich, class = "Enrichment", count =20,  cutoff= 0.05, nchars=70)
dev.off()

# for GC-primed and altered genes
zh_LD_list= unique(filter(zh.all, zfin_id_symbol %in% c(deg.list$`DEG 1`,deg.list$`DEG 2`,deg.list$`DEG 3`))[,2])

try(res_enrich <-disease_enrichment( entities =zh_LD_list, vocabulary = "HGNC",
                                     database = "CURATED" ) # "DEG1+2+3"
    ,silent =T)
table.dis.LD <- res_enrich@qresult
table.dis.LD.sig = table.dis.LD %>% arrange(desc(gg)) %>% filter(FDR <0.01) %>% filter(Count >= 20)

a= str_split(table.dis.LD.sig[,4],pattern = "/")
if (length(a)>1) {
  b= data.frame(Reduce(rbind,a))
} else {b= data.frame(a[[1]][1],a[[1]][2])  
}

table.dis.LD.sig[,"n.input"] = as.numeric(b[,2])

a= str_split(table.dis.LD.sig[,5],pattern = "/")
if (length(a)>1) {
  b= data.frame(Reduce(rbind,a))
} else {b= data.frame(a[[1]][1],a[[1]][2])  
}
table.dis.LD.sig[,"gene.in.catg"] = as.numeric(b[,1])

table.dis.LD.sig[,"ratio.catg"] = table.dis.LD.sig[,"Count"]/table.dis.LD.sig[,"gene.in.catg"]

write.csv(table.dis.LD.sig,"./outputs/s.table6_DisGeNET_enrichment_DEG123.csv")

tiff("./figures/Fig6A_top.disgenet_DEG123.tiff", width = 13, height = 12,units = "cm", res = 1200, compression = "lzw")
plot(res_enrich, class = "Enrichment", count =40,  cutoff= 0.01, nchars=70)
dev.off()

#Venn diagram
dis.name= unique(table.dis.LD.sig$Description)

dep.genes = unique(c(unlist(str_split(table.dis.LD.sig$shared_symbol[c(which(table.dis.LD.sig$Description %in% dis.name[c(4,5)]))],pattern = ";"))) )
schiz.genes = unique(c(unlist(str_split(table.dis.LD.sig$shared_symbol[c(which(table.dis.LD.sig$Description %in% dis.name[c(1)]))],pattern = ";"))) ) 
auti.genes = unique(c(unlist(str_split(table.dis.LD.sig$shared_symbol[c(which(table.dis.LD.sig$Description %in% dis.name[c(3)]))],pattern = ";"))) ) 
bipolar.genes = unique(c(unlist(str_split(table.dis.LD.sig$shared_symbol[c(which(table.dis.LD.sig$Description %in% dis.name[c(2)]))],pattern = ";"))) ) 

dep.disig.zh = zh.all %>% filter(human.Symbol %in% dep.genes)
schiz.disig.zh = zh.all %>% filter(human.Symbol %in% schiz.genes)
auti.disig.zh = zh.all %>% filter(human.Symbol %in% auti.genes)
bipolar.disig.zh = zh.all %>% filter(human.Symbol %in% bipolar.genes)

dep.disig.LD = deg.cpm8 %>% filter(zfin_id_symbol %in% intersect(ab.LD.res.deg.gene,unique(dep.disig.zh$zfin_id_symbol))) 
schiz.disig.LD = deg.cpm8 %>% filter(zfin_id_symbol %in% intersect(ab.LD.res.deg.gene,unique(schiz.disig.zh$zfin_id_symbol)))
auti.disig.LD = deg.cpm8 %>% filter(zfin_id_symbol %in% intersect(ab.LD.res.deg.gene,unique(auti.disig.zh$zfin_id_symbol)))
bipolar.disig.LD = deg.cpm8 %>% filter(zfin_id_symbol %in% intersect(ab.LD.res.deg.gene,unique(bipolar.disig.zh$zfin_id_symbol)))

###

disig.list = list(
  "depressive" = unique(dep.disig.LD$zfin_id_symbol),
  "schizo." = unique(schiz.disig.LD$zfin_id_symbol),
  "Autistic" = unique(auti.disig.LD$zfin_id_symbol),
  "bipolar" = unique(bipolar.disig.LD$zfin_id_symbol)
)


library(eulerr)


fit4= euler(disig.list,shape = "ellipse") # 10 genes between dep. and schizo. are missing because of the adjustment. manually adjust.
venncol4 = viridis(4)

disnet.all.venn =plot(fit4, 
                      quantities = T,
                      fill = venncol4, alpha =0.3,
                      lty = 1,
                      labels = list(font = 4)
)

'
                                    original fitted residuals regionError
depressive                                18     18         0       0.002
schizo.                                  100    100         0       0.013
Autistic                                  39     39         0       0.005
bipolar                                   36     36         0       0.005
depressive&schizo.                        10      0        10       0.035
depressive&Autistic                        4      4         0       0.001
depressive&bipolar                        10     10         0       0.001
schizo.&Autistic                           8      8         0       0.001
schizo.&bipolar                           28     28         0       0.004
Autistic&bipolar                           1      1         0       0.000
depressive&schizo.&Autistic                0      0         0       0.000
depressive&schizo.&bipolar                14     14         0       0.002
depressive&Autistic&bipolar                3      3         0       0.000
schizo.&Autistic&bipolar                   5      5         0       0.001
depressive&schizo.&Autistic&bipolar        9      9         0       0.001

diagError: 0.035 
stress:    0.007
'

# need to add losses in fitted  depressive&schizo-> 10
tiff("./figures/Fig6B_disgnet.all.venn4.12x12.tiff", height = 12,width = 12, res=1200, units = "cm", compression = "lzw")
plot(disnet.all.venn)
dev.off()


####log2fc
for (i in c(5:8)) {
  deg = get(paste0("deg.cpm",i))
  dep.disig = deg %>% filter(zfin_id_symbol %in% unique(dep.disig.LD$zfin_id_symbol)) 
  schiz.disig = deg %>% filter(zfin_id_symbol %in% unique(schiz.disig.LD$zfin_id_symbol))
  auti.disig = deg %>% filter(zfin_id_symbol %in% unique(auti.disig.LD$zfin_id_symbol))
  bipolar.disig= deg %>% filter(zfin_id_symbol %in% unique(bipolar.disig.LD$zfin_id_symbol))
  name1=paste0("dep.disig",i)
  name2=paste0("schiz.disig",i)
  name3=paste0("auti.disig",i)
  name4=paste0("bipolar.disig",i)
  assign(name1,dep.disig)
  assign(name2,schiz.disig)
  assign(name3,auti.disig)
  assign(name4,bipolar.disig)
  
}


#### intersection between LD-deg and disgenet enriched DEG
dp.int.ld = intersect(disig.list[["depressive"]],ab.LD.res.deg.gene)
length(dp.int.ld)

sch.int.ld = intersect(disig.list[["schizo."]],ab.LD.res.deg.gene)
length(sch.int.ld)

aut.int.ld = intersect(disig.list[["Autistic"]],ab.LD.res.deg.gene) 
length(aut.int.ld)

bp.int.ld = intersect(disig.list[["bipolar"]],ab.LD.res.deg.gene) 
length(bp.int.ld)

LD_DEG_disg= data.frame("LD_DEG"= c(dp.int.ld,sch.int.ld,aut.int.ld,bp.int.ld), "source"= c(rep("depressive",length(dp.int.ld)),
                                                                                        rep("schizo.",length(sch.int.ld)),
                                                                                        rep("Autistic",length(aut.int.ld)),
                                                                                        rep("bipolar",length(bp.int.ld))
                                                                                        )
                        )

write.csv(LD_DEG_disg, "./outputs/s.table7_disease.associated.LD-DEG.csv")

####log2fc
ht.mt = data.frame(dep.disig5$logFC,dep.disig6$logFC,
                   dep.disig7$logFC,dep.disig8$logFC)

row.names(ht.mt)=dep.disig5$zfin_id_symbol

ht.FDR= data.frame(dep.disig5$FDR,dep.disig6$FDR,
                   dep.disig7$FDR,dep.disig8$FDR)

row.names(ht.FDR)=dep.disig5$zfin_id_symbol

ht.FDR[ht.FDR <= 0.05] <-"*"
ht.FDR[ht.FDR > 0.05] <-NA


library (gplots)
breaks=seq(-3,3,0.01)
ycol = colorRampPalette(c("#313695","#4575B4","#74ADD1","#ABD9E9" ,"#E0F3F8","#FFFFBF" ,"#FEE090" ,"#FDAE61", 
                           "#F46D43","#D73027","#A50026" ))(n=600)


st_sig_tx_1= as.matrix(ht.mt)
distfunc <- function(x) dist(x, method="euclidean")
hclustfunc <- function(x) hclust(x, method="ward.D2")
dep.rowname = row.names(ht.mt)
dep.rowname[dep.rowname %in% dp.int.ld] <- "red"
dep.rowname[dep.rowname != "red"] <- "black"
ldrowcols=dep.rowname

tiff(paste0("./figures/Fig6d_dp.heatmap.10x5_a",nrow(st_sig_tx_1) ,".tiff"), height = 10,width = 5, res=1200, units = "cm", compression = "lzw")

try(heatmap.2((st_sig_tx_1), col=mycol, scale='none',
               dendrogram = "none",Colv = F,breaks=breaks,hclustfun=hclustfunc, distfun=distfunc,
               density='density', trace = 'none',
               na.color = F,
               offsetRow = 0,
               offsetCol = c(0,0.2),
               cexCol = 1,
               colsep = 4,
               labCol = c("6","13","120","LD"),
               labRow = "",
               keysize = 1 ), silent = T)

dev.off()

ht.mt = data.frame(schiz.disig5$logFC,schiz.disig6$logFC,
                   schiz.disig7$logFC,schiz.disig8$logFC)

row.names(ht.mt)=schiz.disig5$zfin_id_symbol

ht.FDR= data.frame(schiz.disig5$FDR,schiz.disig6$FDR,
                   schiz.disig7$FDR,schiz.disig8$FDR)

row.names(ht.FDR)=schiz.disig5$zfin_id_symbol

ht.FDR[ht.FDR <= 0.05] <-"*"
ht.FDR[ht.FDR > 0.05] <-NA


library (gplots)
breaks=seq(-3,3,0.01)
mycol = colorRampPalette(c("#313695","#4575B4","#74ADD1","#ABD9E9" ,"#E0F3F8","#FFFFBF" ,"#FEE090" ,"#FDAE61", 
                           "#F46D43","#D73027","#A50026" ))(n=600)


st_sig_tx_1= as.matrix(ht.mt)



distfunc <- function(x) dist(x, method="euclidean")
hclustfunc <- function(x) hclust(x, method="ward.D2")

tiff(paste0("./figures/Fig6d_sch.heatmap.10x5_",nrow(st_sig_tx_1),".tiff"), height = 10,width = 5, res=1200, units = "cm", compression = "lzw")

try(heatmap.2((st_sig_tx_1), col=mycol, scale='none',
               dendrogram = "none",Colv = F,breaks=breaks,hclustfun=hclustfunc, distfun=distfunc,
               density='density', trace = 'none',
               na.color = F,
               offsetRow = 0,
               offsetCol = c(0,0.2),
               cexCol = 1,
               colsep = 4,
               labCol = c("6","13","120","LD"),
               labRow = "",
               keysize = 1 ),silent = T)

dev.off()


ht.mt = data.frame(auti.disig5$logFC,auti.disig6$logFC,
                   auti.disig7$logFC,auti.disig8$logFC)

row.names(ht.mt)=auti.disig5$zfin_id_symbol

ht.FDR= data.frame(auti.disig5$FDR,auti.disig6$FDR,
                   auti.disig7$FDR,auti.disig8$FDR)

row.names(ht.FDR)=auti.disig5$zfin_id_symbol

ht.FDR[ht.FDR <= 0.05] <-"*"
ht.FDR[ht.FDR > 0.05] <-NA


library (gplots)
breaks=seq(-3,3,0.01)
mycol = colorRampPalette(c("#313695","#4575B4","#74ADD1","#ABD9E9" ,"#E0F3F8","#FFFFBF" ,"#FEE090" ,"#FDAE61", 
                           "#F46D43","#D73027","#A50026" ))(n=600)


st_sig_tx_1= as.matrix(ht.mt)

distfunc <- function(x) dist(x, method="euclidean")
hclustfunc <- function(x) hclust(x, method="ward.D2")

tiff(paste0("./figures/Fig6d_auti.heatmap.10x5_",nrow(st_sig_tx_1),".tiff"), height = 10,width = 5, res=1200, units = "cm", compression = "lzw")

try(heatmap.2((st_sig_tx_1), col=mycol, scale='none',
               dendrogram = "none",Colv = F,breaks=breaks,hclustfun=hclustfunc, distfun=distfunc,
               #main=c("oxt_regulation, log2FC"),
               density='density', trace = 'none',
               na.color = F,
               offsetRow = 0,
               offsetCol = c(0,0.2),
               cexCol = 1,
               colsep = 4,
               labCol = c("6","13","120","LD"),
               labRow = "",
               keysize = 1 ), silent = T)

dev.off()

ht.mt = data.frame(bipolar.disig5$logFC,bipolar.disig6$logFC,
                   bipolar.disig7$logFC,bipolar.disig8$logFC)

row.names(ht.mt)=bipolar.disig5$zfin_id_symbol

ht.FDR= data.frame(bipolar.disig5$FDR,bipolar.disig6$FDR,
                   bipolar.disig7$FDR,bipolar.disig8$FDR)

row.names(ht.FDR)=bipolar.disig5$zfin_id_symbol

ht.FDR[ht.FDR <= 0.05] <-"*"
ht.FDR[ht.FDR > 0.05] <-NA


library (gplots)
breaks=seq(-3,3,0.01)
mycol = colorRampPalette(c("#313695","#4575B4","#74ADD1","#ABD9E9" ,"#E0F3F8","#FFFFBF" ,"#FEE090" ,"#FDAE61", 
                           "#F46D43","#D73027","#A50026" ))(n=600)


st_sig_tx_1= as.matrix(ht.mt)



distfunc <- function(x) dist(x, method="euclidean")
hclustfunc <- function(x) hclust(x, method="ward.D2")

tiff(paste0("./figures/Fig6d_bipolar.heatmap.10x5_",nrow(st_sig_tx_1),".tiff"), height = 10,width = 5, res=1200, units = "cm", compression = "lzw")

try(heatmap.2((st_sig_tx_1), col=mycol, scale='none',
               dendrogram = "none",Colv = F,breaks=breaks,hclustfun=hclustfunc, distfun=distfunc,
               density='density', trace = 'none',
               na.color = F,
               offsetRow = 0,
               offsetCol = c(0,0.2),
               cexCol = 1,
               colsep = 4,
               labCol = c("6","13","120","LD"),
               labRow = "",
               keysize = 1 ),silent = T)

dev.off()


#### all disegenet


top_disLD= unique(c(schiz.disig$ensembl_gene_id,dep.disig$ensembl_gene_id, bipolar.disig$ensembl_gene_id,auti.disig$ensembl_gene_id))

for (i in c(5:8)) {
  deg = get(paste0("deg.cpm",i))
  dep.disig = deg %>% filter(ensembl_gene_id %in% top_disLD) 
  name1=paste0("top.disig",i)
  assign(name1,dep.disig)
 
}

ht.mt = data.frame(top.disig5$logFC,top.disig6$logFC,
                   top.disig7$logFC,top.disig8$logFC)

row.names(ht.mt)=top.disig5$zfin_id_symbol

ht.FDR= data.frame(top.disig5$FDR,top.disig6$FDR,
                   top.disig7$FDR,top.disig8$FDR)

row.names(ht.FDR)=top.disig5$zfin_id_symbol

ht.FDR[ht.FDR <= 0.05] <-"*"
ht.FDR[ht.FDR > 0.05] <-NA

ht.mt[ht.FDR == "NA"] <-NA

library (gplots)
breaks=seq(-3,3,0.01)
mycol = colorRampPalette(c("#313695","#4575B4","#74ADD1","#ABD9E9" ,"#E0F3F8","#FFFFBF" ,"#FEE090" ,"#FDAE61", 
                           "#F46D43","#D73027","#A50026" ))(n=600)


st_sig_tx_1= as.matrix(ht.mt)


distfunc <- function(x) dist(x, method="euclidean")
hclustfunc <- function(x) hclust(x, method="ward.D2")

hc = hclustfunc(distfunc(st_sig_tx_1))
ht.label.order = hc$labels[hc$order]
write(ht.label.order,"./figures/Fig6e.top.dis.LD.DEGS.label.order.txt")

rhc = hclust(dist(st_sig_tx_1, method="euclidean"), method="ward.D2")
gr.row <- cutree(rhc, 4) 

desturate =viridis(4)

tiff(paste0("./figures/Fig6e_top.dis.LD.DEGS.heatmap.a",nrow(st_sig_tx_1) ,".tiff"), height = 10,width = 5, res=1200, units = "cm", compression = "lzw")

try(heatmap.2((st_sig_tx_1), col=mycol, scale='none',
              dendrogram = "row",Colv = F,breaks=breaks,hclustfun=hclustfunc, distfun=distfunc,
              density='density', trace = 'none',
              na.color = T,
              offsetRow = 0,
              offsetCol = 0,
              cexCol = 1,
              cexRow = 1,
              #colsep = 4,
              labCol = c("6","13","120","LD"),
              RowSideColors=desturate[gr.row],
              lwid = c(1,2.5),
              
              #cellnote=adult_specific.ht.FDR,
              #notecex=1,
              #notecol="black",
              labRow = "",
              keysize = 0.5 ), silent = T)

dev.off()
#### expression trajectory, consistent or de novo

#selection
for (i in c(5:8)) {
  deg = get(paste0("deg.cpm",i))
  deg = deg %>% filter(zfin_id_symbol %in% unique(LD_DEG_disg$LD_DEG))
  name=paste0("LD.all.disig",i)
  assign(name,deg)
}


#trans_temporal
main.LD.disig.FDR = Reduce(intersect, list(LD.all.disig5$zfin_id_symbol[c(which(LD.all.disig5$FDR <0.01))],
                                       LD.all.disig6$zfin_id_symbol[c(which(LD.all.disig6$FDR <0.01))],
                                       LD.all.disig7$zfin_id_symbol[c(which(LD.all.disig7$FDR <0.01))],
                                       LD.all.disig8$zfin_id_symbol[c(which(LD.all.disig8$FDR <0.01))]))

matched.alt = c()
  
  for (i in 1:nrow(LD.all.disig8)) {
    if (abs(LD.all.disig8$logFC[i])> log2(2) &
        abs(LD.all.disig7$logFC[i])> log2(1.5) &
        (abs(LD.all.disig6$logFC[i])> log2(2) |
        abs(LD.all.disig5$logFC[i])> log2(2))){
      if (LD.all.disig8$logFC[i]*LD.all.disig7$logFC[i] > 0 & 
         LD.all.disig8$logFC[i]*LD.all.disig6$logFC[i] > 0 &
         LD.all.disig8$logFC[i]*LD.all.disig5$logFC[i] > 0){
       matched.alt[i]=1
      } else if (LD.all.disig8$logFC[i]*LD.all.disig7$logFC[i] > 0 & 
                 (LD.all.disig8$logFC[i]*LD.all.disig6$logFC[i] > 0 |
                 LD.all.disig8$logFC[i]*LD.all.disig5$logFC[i] > 0)){
        matched.alt[i]=1
      } else  {
       matched.alt[i]=0
      }
    }else { 
      matched.alt[i]=0
      }
}
 

trans_temporal.alt= LD.all.disig8[c(which(matched.alt ==1)),8]
trans_temporal.alt.sig = intersect(trans_temporal.alt,main.LD.disig.FDR)


adult_specific.alt = c()

for (i in 1:nrow(LD.all.disig8)){
  if (abs(LD.all.disig8$logFC[i])> log2(2) &
      abs(LD.all.disig6$logFC[i])< log2(1.5) &
       abs(LD.all.disig5$logFC[i])< log2(1.5)){
    adult_specific.alt[i]=1
  }else { 
    adult_specific.alt[i]=0
  }
}

adult_specific.alt= LD.all.disig8[c(which(adult_specific.alt == 1)),8]
adult_specific.alt.sig = intersect(adult_specific.alt,LD.all.disig8$zfin_id_symbol[c(which(LD.all.disig8$FDR <0.01))])

#human_orthologs
combi.dis.LD_deg = c(trans_temporal.alt.sig,adult_specific.alt.sig)
combi.dis.LD_deg.zh = filter(zh.all,zfin_id_symbol %in% combi.dis.LD_deg)
#write.csv(combi.dis.LD_deg.zh,"./outputs/combi.dis.LD_deg.zh.csv")

####heatmap

maint.ht.mt = cbind("d6"=LD.all.disig5$logFC[c(which(LD.all.disig5$zfin_id_symbol %in% trans_temporal.alt.sig))],
                    "d13"=LD.all.disig6$logFC[c(which(LD.all.disig6$zfin_id_symbol %in% trans_temporal.alt.sig))],
                    "d120"=LD.all.disig7$logFC[c(which(LD.all.disig7$zfin_id_symbol %in% trans_temporal.alt.sig))],
                    "LD"=LD.all.disig8$logFC[c(which(LD.all.disig8$zfin_id_symbol %in% trans_temporal.alt.sig))])

st_sig_tx_1= as.matrix(maint.ht.mt)
distfunc <- function(x) dist(x, method="euclidean")
hclustfunc <- function(x) hclust(x, method="ward.D2")
row.names(st_sig_tx_1)= LD.all.disig8$zfin_id_symbol[c(which(LD.all.disig5$zfin_id_symbol %in% trans_temporal.alt.sig))]

hc = hclustfunc(distfunc(st_sig_tx_1))
ht.label.order = hc$labels[hc$order]
write(ht.label.order,"./figures/Fig6e_trans_temporal.dis.DEGS.label.order.txt")
nrow(maint.ht.mt)#60
tiff(paste0("./figures/s.Fig6e_trans_temporal.dis.DEGS.heatmap_a",nrow(st_sig_tx_1) ,".tiff"), height = 7,width = 5, res=1200, units = "cm", compression = "lzw")
try(heatmap.2((st_sig_tx_1), col=mycol, scale='none',
              dendrogram = "row",Colv = F,breaks=breaks,hclustfun=hclustfunc, distfun=distfunc,
              density='density', trace = 'none',
              na.color = F,
              offsetRow = 0,
              offsetCol = c(0,0),
              cexCol = 1,
              cexRow = 0.4,
              colsep = 4,
              labCol = c("6","13","120","LD"),
              labRow = "",
              keysize = 1 ), silent = T)

dev.off()


adult_specific.ht.mt = cbind("d6"=LD.all.disig5$logFC[c(which(LD.all.disig5$zfin_id_symbol %in% adult_specific.alt.sig))],
                    "d13"=LD.all.disig6$logFC[c(which(LD.all.disig6$zfin_id_symbol %in% adult_specific.alt.sig))],
                    "d120"=LD.all.disig7$logFC[c(which(LD.all.disig7$zfin_id_symbol %in% adult_specific.alt.sig))],
                    "LD"=LD.all.disig8$logFC[c(which(LD.all.disig8$zfin_id_symbol %in% adult_specific.alt.sig))])


st_sig_tx_1= as.matrix(adult_specific.ht.mt)
distfunc <- function(x) dist(x, method="euclidean")
hclustfunc <- function(x) hclust(x, method="ward.D2")
row.names(st_sig_tx_1)= LD.all.disig8$zfin_id_symbol[c(which(LD.all.disig8$zfin_id_symbol %in% adult_specific.alt.sig))]

hc = hclustfunc(distfunc(st_sig_tx_1))
ht.label.order = hc$labels[hc$order]
write(ht.label.order,"./figures/Fig6e_adult_specific.dis.DEGS.label.order.txt")
nrow(adult_specific.ht.mt)#65


tiff(paste0("./figures/Fig6e_adult_specific.dis.DEGS.heatmap.4x4_",nrow(st_sig_tx_1) ,".tiff"), height = 7,width = 5, res=1200, units = "cm", compression = "lzw")

try(heatmap.2((st_sig_tx_1), col=mycol, scale='none',
              dendrogram = "row",Colv = F,breaks=breaks,hclustfun=hclustfunc, distfun=distfunc,
              density='density', trace = 'none',
              na.color = F,
              offsetRow = 0,
              offsetCol = c(0,0.2),
              cexCol = 1,
              cexRow = 0.8,
              colsep = 4,
              labCol = c("6","13","120","LD"),
              #cellnote=adult_specific.ht.FDR,
              #notecex=1,
              #notecol="black",
              labRow = "",
              keysize = 1 ), silent = T)

dev.off()

#combi
maint.adult_specific.ht.mt = rbind(maint.ht.mt,adult_specific.ht.mt)

st_sig_tx_1= as.matrix(maint.adult_specific.ht.mt)
distfunc <- function(x) dist(x, method="euclidean")
hclustfunc <- function(x) hclust(x, method="ward.D2")
row.names(st_sig_tx_1)= c(LD.all.disig8$zfin_id_symbol[c(which(LD.all.disig5$zfin_id_symbol %in% trans_temporal.alt.sig))],
  LD.all.disig8$zfin_id_symbol[c(which(LD.all.disig8$zfin_id_symbol %in% adult_specific.alt.sig))]
  )

hc = hclustfunc(distfunc(st_sig_tx_1))
ht.label.order = hc$labels[hc$order]
write(ht.label.order,"./figures/Fig6e.combi.dis.DEGS.label.order.txt")

rhc = hclust(dist(st_sig_tx_1, method="euclidean"), method="ward.D2")
gr.row <- cutree(rhc, 2) 

desturate =viridis(2)

tiff(paste0("./figures/Fig6e_dis.combi.dis.DEGS.heatmap.",nrow(st_sig_tx_1) ,".tiff"), height = 10,width = 5, res=1200, units = "cm", compression = "lzw")

try(heatmap.2((st_sig_tx_1), col=mycol, scale='none',
              dendrogram = "row",Colv = F,breaks=breaks,hclustfun=hclustfunc, distfun=distfunc,
              density='density', trace = 'none',
              na.color = F,
              offsetRow = 0,
              offsetCol = 0,
              cexCol = 1,
              cexRow = 1,
              #colsep = 4,
              labCol = c("6","13","120","LD"),
              RowSideColors=desturate[gr.row],
              lwid = c(1,2.5),
              
              #cellnote=adult_specific.ht.FDR,
              #notecex=1,
              #notecol="black",
              labRow = "",
              keysize = 0.5 ), silent = T)

dev.off()

all_dis_LD_deg=rbind(LD.all.disig5,LD.all.disig6,LD.all.disig7,LD.all.disig8)
#write.csv(all_dis_LD_deg,"./outputs/all_dis_LD_deg.csv")
#write.csv(LD.all.disig7,"./outputs/LD.all.disig7.csv")
#write.csv(LD.all.disig8,"./outputs/LD.all.disig8.csv")

##### version for ALL DEGs


#selection
for (i in c(5:8)) {
  deg = get(paste0("deg.cpm",i))
  deg = deg %>% filter(zfin_id_symbol %in% unique(c(deg.list$`DEG 1`,deg.list$`DEG 2`,deg.list$`DEG 3`)))
  name=paste0("LD.all.",i)
  assign(name,deg)
}
deg12.primed

#trans_temporal
main.LD.all.FDR = Reduce(intersect, list(LD.all.5$zfin_id_symbol[c(which(LD.all.5$FDR <0.01))],
                                           LD.all.6$zfin_id_symbol[c(which(LD.all.6$FDR <0.01))],
                                           LD.all.7$zfin_id_symbol[c(which(LD.all.7$FDR <0.05))],
                                           LD.all.8$zfin_id_symbol[c(which(LD.all.8$FDR <0.01))]))

matched.alt = c()

for (i in 1:nrow(LD.all.8)) {
  if (abs(LD.all.8$logFC[i])> log2(2) &
      abs(LD.all.7$logFC[i])> log2(1.5) &
      (abs(LD.all.6$logFC[i])> log2(2) |
       abs(LD.all.5$logFC[i])> log2(2))){
    if (LD.all.8$logFC[i]*LD.all.7$logFC[i] > 0 & 
        LD.all.8$logFC[i]*LD.all.6$logFC[i] > 0 &
        LD.all.8$logFC[i]*LD.all.5$logFC[i] > 0){
      matched.alt[i]=1
    } else if (LD.all.8$logFC[i]*LD.all.7$logFC[i] > 0 & 
               (LD.all.8$logFC[i]*LD.all.6$logFC[i] > 0 |
                LD.all.8$logFC[i]*LD.all.5$logFC[i] > 0)){
      matched.alt[i]=1
    } else  {
      matched.alt[i]=0
    }
  }else { 
    matched.alt[i]=0
  }
}


trans_temporal.alt= LD.all.8[c(which(matched.alt ==1)),8]
trans_temporal.alt.sig = intersect(trans_temporal.alt,main.LD.all.FDR)
trans_temporal.alt.sig.table = dplyr::filter(LD.all.8, zfin_id_symbol %in% trans_temporal.alt.sig )
write.csv(trans_temporal.alt.sig.table, "./outputs/s.table8.trans_temporal_all.LD_DEGs.csv")

adult_specific.alt = c()

for (i in 1:nrow(LD.all.8)){
  if (abs(LD.all.8$logFC[i])> log2(2) & 
      abs(LD.all.6$logFC[i])< log2(1.5) &
      abs(LD.all.5$logFC[i])< log2(1.5)){
    adult_specific.alt[i]=1
  }else { 
    adult_specific.alt[i]=0
  }
}

adult_specific.alt= LD.all.8[c(which(adult_specific.alt == 1)),8]
adult_specific.alt.sig = intersect(adult_specific.alt,LD.all.8$zfin_id_symbol[c(which(LD.all.8$FDR <0.01))])
adult_specific.alt.sig.table = dplyr::filter(LD.all.8, zfin_id_symbol %in% adult_specific.alt.sig )
write.csv(adult_specific.alt.sig.table,"./outputs/s.table8.adult_specific_all.LD_DEGs.csv")

#human_orthologs
combi.LD_deg = c(trans_temporal.alt.sig,adult_specific.alt.sig)
combi.LD_deg.zh = dplyr::filter(zh.all, zfin_id_symbol %in% combi.LD_deg)
#write.csv(combi.dis.LD_deg.zh,"./outputs/combi.dis.LD_deg.zh.csv")

####heatmap trans_temporal.alt.sig

maint.ht.mt = cbind("d6"=LD.all.5$logFC[c(which(LD.all.5$zfin_id_symbol %in% trans_temporal.alt.sig))],
                    "d13"=LD.all.6$logFC[c(which(LD.all.6$zfin_id_symbol %in% trans_temporal.alt.sig))],
                    "d120"=LD.all.7$logFC[c(which(LD.all.7$zfin_id_symbol %in% trans_temporal.alt.sig))],
                    "LD"=LD.all.8$logFC[c(which(LD.all.8$zfin_id_symbol %in% trans_temporal.alt.sig))])

st_sig_tx_1= as.matrix(maint.ht.mt)
distfunc <- function(x) dist(x, method="euclidean")
hclustfunc <- function(x) hclust(x, method="ward.D2")
row.names(st_sig_tx_1)= LD.all.8$zfin_id_symbol[c(which(LD.all.5$zfin_id_symbol %in% trans_temporal.alt.sig))]

hc = hclustfunc(distfunc(st_sig_tx_1))
ht.label.order = hc$labels[hc$order]
write(ht.label.order,"./figures/Fig6e_trans_temporal.LD.DEGS2.label.order.txt")
nrow(maint.ht.mt)#560
tiff(paste0("./figures/Fig6e_trans_temporal.LD.DEGS2.heatmap_a",nrow(st_sig_tx_1) ,".tiff"), height = 7,width = 5, res=1200, units = "cm", compression = "lzw")
try(heatmap.2((st_sig_tx_1), col=mycol, scale='none',
              dendrogram = "row",Colv = F,breaks=breaks,hclustfun=hclustfunc, distfun=distfunc,
              density='density', trace = 'none',
              na.color = F,
              offsetRow = 0,
              offsetCol = c(0,0),
              cexCol = 1,
              cexRow = 0.4,
              colsep = 4,
              labCol = c("6","13","120","LD"),
              labRow = "",
              keysize = 1 ), silent = T)

dev.off()

#heatmap for adult_specific.alt.sig

adult_specific.ht.mt = cbind("d6"=LD.all.5$logFC[c(which(LD.all.5$zfin_id_symbol %in% adult_specific.alt.sig))],
                     "d13"=LD.all.6$logFC[c(which(LD.all.6$zfin_id_symbol %in% adult_specific.alt.sig))],
                     "d120"=LD.all.7$logFC[c(which(LD.all.7$zfin_id_symbol %in% adult_specific.alt.sig))],
                     "LD"=LD.all.8$logFC[c(which(LD.all.8$zfin_id_symbol %in% adult_specific.alt.sig))])


st_sig_tx_1= as.matrix(adult_specific.ht.mt)
distfunc <- function(x) dist(x, method="euclidean")
hclustfunc <- function(x) hclust(x, method="ward.D2")
row.names(st_sig_tx_1)= LD.all.8$zfin_id_symbol[c(which(LD.all.8$zfin_id_symbol %in% adult_specific.alt.sig))]

hc = hclustfunc(distfunc(st_sig_tx_1))
ht.label.order = hc$labels[hc$order]
write(ht.label.order,"./figures/Fig6e_adult_specific.ELSprimed.LD.DEGS.label.order.txt")
nrow(adult_specific.ht.mt)#69


tiff(paste0("./figures/Fig6e_adult_specific.ELSprimed.LD.DEGS.2.heatmap.4x4_a",nrow(st_sig_tx_1) ,".tiff"),
     height = 7,width = 5, res=1200, units = "cm", compression = "lzw")

try(heatmap.2((st_sig_tx_1), col=mycol, scale='none',
              dendrogram = "row",Colv = F,breaks=breaks,hclustfun=hclustfunc, distfun=distfunc,
              density='density', trace = 'none',
              na.color = F,
              offsetRow = 0,
              offsetCol = c(0,0.2),
              cexCol = 1,
              cexRow = 0.8,
              colsep = 4,
              labCol = c("6","13","120","LD"),
              #cellnote=adult_specific.ht.FDR,
              #notecex=1,
              #notecol="black",
              labRow = "",
              keysize = 1 ), silent = T)

dev.off()

#combi
maint.adult_specific.ht.mt = rbind(maint.ht.mt,adult_specific.ht.mt)

st_sig_tx_1= as.matrix(maint.adult_specific.ht.mt)
distfunc <- function(x) dist(x, method="euclidean")
hclustfunc <- function(x) hclust(x, method="ward.D2")
row.names(st_sig_tx_1)= c(LD.all.8$zfin_id_symbol[c(which(LD.all.5$zfin_id_symbol %in% trans_temporal.alt.sig))],
                          LD.all.8$zfin_id_symbol[c(which(LD.all.8$zfin_id_symbol %in% adult_specific.alt.sig))]
)

hc = hclustfunc(distfunc(st_sig_tx_1))
ht.label.order = hc$labels[hc$order]
write(ht.label.order,"./figures/Fig6e.ELSprimed.LD.DEGS.2.label.order.txt")

rhc = hclust(dist(st_sig_tx_1, method="euclidean"), method="ward.D2")
gr.row <- cutree(rhc, 4)
tree.table = data.frame("gene" = names(gr.row), "cluster"= gr.row)
write.csv(tree.table,"./figures/Fig6e.ELSprimed.LD.DEGS.2.cluster.order.csv")

desturate =viridis(4)

tiff(paste0("./figures/Fig6e_ELSprimed.LD.DEGS.2.heatmap.a",nrow(st_sig_tx_1) ,".tiff"), height = 10,width = 5, res=1200, units = "cm", compression = "lzw")

try(heatmap.2((st_sig_tx_1), col=mycol, scale='none',
              dendrogram = "row",Colv = F,breaks=breaks,hclustfun=hclustfunc, distfun=distfunc,
              density='density', trace = 'none',
              na.color = F,
              offsetRow = 0,
              offsetCol = 0,
              cexCol = 1,
              cexRow = 1,
              #colsep = 4,
              labCol = c("6","13","120","LD"),
              #RowSideColors=desturate[gr.row],
              lwid = c(1,2.5),
              
              #cellnote=adult_specific.ht.FDR,
              #notecex=1,
              #notecol="black",
              labRow = "",
              keysize = 0.5 ), silent = T)

dev.off()

all_LD_deg=rbind(LD.all.5,LD.all.6,LD.all.7,LD.all.8)
#write.csv(all_dis_LD_deg,"./outputs/all_dis_LD_deg.csv")
#write.csv(LD.all.7,"./outputs/LD.all.disig7.csv")
#write.csv(LD.all.8,"./outputs/LD.all.disig8.csv")

trans_temporal_sch_dep_aut = tree.table %>% 
  filter(cluster == 2|3) %>% filter(gene %in% c(LD_DEG_disg$LD_DEG[c(which(LD_DEG_disg$source %in% c("depressive", "schizo.", "Autistic")))]))
trans_temporal_sch_dep_exp.d6 = filter(deg.cpm5,zfin_id_symbol %in% trans_temporal_sch_dep_aut$gene )
#write.csv(trans_temporal_sch_dep_exp.d6,"./outputs/trans_temporal_sch_dep_exp.d6.csv")


###cell type

# load enrichr
library(reticulate)
library(plyr)
os <- import("os")
os$listdir(".")

json = import("json")
requests= import("requests")

#get human orthologues
tt_ld= na.omit(zh.all %>% filter(zfin_id_symbol %in%trans_temporal.alt.sig ))
tt_ld=unique(tt_ld[,2])

#transtemporal
ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
genes_str = as.character(paste0(tt_ld,collapse = "\n"))

description = 'transtemporal_DEGs'

payload = list('list'= c("None",genes_str),
            'description'= c("None",description))

response = requests$post(ENRICHR_URL, files=payload)

data = json$loads(response$text)
print(data)

# get enrichment
ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
query_string = '?userListId=%s&backgroundType=%s'
user_list_id = data$userListId
gene_set_library = 'Tabula_Muris'
url = paste0(ENRICHR_URL,"?","userListId=",user_list_id,
             "&backgroundType=",gene_set_library)

response = requests$get(url)

data2 = json$loads(response$text)
print(data2)

### get_result

#check the ID and modifiy the .py file.
source_python("./code/enrichr_get_transtemporal_result.py")

# adult_specific

as_ld= na.omit(zh.all %>% filter(zfin_id_symbol %in%adult_specific.alt.sig ))
as_ld=unique(as_ld[,2])


ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
genes_str = as.character(paste0(as_ld,collapse = "\n"))

description = 'adult_specific_DEGs'

payload = list('list'= c("None",genes_str),
               'description'= c("None",description))

response = requests$post(ENRICHR_URL, files=payload)

data = json$loads(response$text)
print(data)

# get enrichment
ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
query_string = '?userListId=%s&backgroundType=%s'
user_list_id = data$userListId
gene_set_library = 'Tabula_Muris'
url = paste0(ENRICHR_URL,"?","userListId=",user_list_id,
             "&backgroundType=",gene_set_library)

response = requests$get(url)

data2 = json$loads(response$text)
print(data2)

### get_result

source_python("./code/enrichr_get_adult_specific_result.py")


#read_table
tp_celltype_enrich = read.table("./outputs/transtemporal_DEGs123.txt",sep = "\t",header = T)
as_celltype_enrich = read.table("./outputs/adult_specific_DEGs123.txt",sep = "\t",header = T)

sig.tp.cell= filter(tp_celltype_enrich,Adjusted.P.value <0.05)[,c(1,2,3,4,7:9)]
sig.as.cell= filter(as_celltype_enrich,Adjusted.P.value <0.05)[,c(1,2,3,4,7:9)]

sig.tp.cell[,"source"] ="transtemporal"
sig.as.cell[,"source"] ="adult-specific"

sig.cell_type123 = rbind(sig.tp.cell,sig.as.cell)

write.csv(sig.cell_type123,"./outputs/s.table9.sig.cell_type.csv")

###for DEG1+2

#get human orthologues
tt_ld= na.omit(zh.all %>% filter(zfin_id_symbol %in%trans_temporal.alt.sig ) %>% filter(zfin_id_symbol%in% deg12.primed))
tt_ld=unique(tt_ld[,2])

#transtemporal
ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
genes_str = as.character(paste0(tt_ld,collapse = "\n"))

description = 'transtemporal_DEGs'

payload = list('list'= c("None",genes_str),
               'description'= c("None",description))

response = requests$post(ENRICHR_URL, files=payload)

data = json$loads(response$text)
print(data)

# get enrichment
ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
query_string = '?userListId=%s&backgroundType=%s'
user_list_id = data$userListId
gene_set_library = 'Tabula_Muris'
url = paste0(ENRICHR_URL,"?","userListId=",user_list_id,
             "&backgroundType=",gene_set_library)

response = requests$get(url)

data2 = json$loads(response$text)
print(data2)

### get_result

#check the ID and modifiy the .py file.
source_python("./code/enrichr_get_transtemporal_result.py")

# adult_specific

as_ld= na.omit(zh.all %>% filter(zfin_id_symbol %in%adult_specific.alt.sig )%>% filter(zfin_id_symbol%in% deg12.primed))
as_ld=unique(as_ld[,2])


ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
genes_str = as.character(paste0(as_ld,collapse = "\n"))

description = 'adult_specific_DEGs'

payload = list('list'= c("None",genes_str),
               'description'= c("None",description))

response = requests$post(ENRICHR_URL, files=payload)

data = json$loads(response$text)
print(data)

# get enrichment
ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
query_string = '?userListId=%s&backgroundType=%s'
user_list_id = data$userListId
gene_set_library = 'Tabula_Muris'
url = paste0(ENRICHR_URL,"?","userListId=",user_list_id,
             "&backgroundType=",gene_set_library)

response = requests$get(url)

data2 = json$loads(response$text)
print(data2)

### get_result

source_python("./code/enrichr_get_adult_specific_result.py")


#read_table
tp_celltype_enrich = read.table("./outputs/transtemporal_DEGs.txt",sep = "\t",header = T)
as_celltype_enrich = read.table("./outputs/adult_specific_DEGs.txt",sep = "\t",header = T)

sig.tp.cell= filter(tp_celltype_enrich,Adjusted.P.value <0.05)[,c(1,2,3,4,7:9)]
sig.as.cell= filter(as_celltype_enrich,Adjusted.P.value <0.05)[,c(1,2,3,4,7:9)]

sig.tp.cell[,"source"] ="transtemporal"
sig.as.cell[,"source"] ="adult-specific"

sig.cell_type12 = rbind(sig.tp.cell,sig.as.cell)

write.csv(sig.cell_type12,"./outputs/s.table9.sig.cell_type12.csv")

# epigenetic regulators ---------------------------------------------------

#downloaded tables from QuickGO with keywords,"epigenetic regulation of gene expression", "DNA modification", "RNA modification","RNA processing", "histone modification"

epi_reg = unique(read.csv("./data/QuickGO_annotations_epigenetic_reg.tsv", sep="\t",header = T)[,3])
DNA_modi = unique(read.csv("./data/QuickGO_annotations_DNA_modification.tsv", sep="\t",header = T)[,3])
hist_modi = unique(read.csv("./data/QuickGO_annotations_histone_modification.tsv", sep="\t",header = T)[,3])
RNA_proc = unique(read.csv("./data/QuickGO_annotations_RNA_processing.tsv", sep="\t",header = T)[,3])
RNA_modi = unique(read.csv("./data/QuickGO_annotations_RNA_modification.tsv", sep="\t",header = T)[,3])


#selection

epi_terms = c("epi_reg","DNA_modi","hist_modi", "RNA_modi","RNA_proc")

for (i in c(5:8)) {
  deg = get(paste0("LD.all.",i))
  deg = deg %>% filter(mean.cpm.pos > 5|mean.cpm.wt >5)%>% filter(t.test < 0.05)%>%filter(FDR < 0.01)%>%
    filter(abs(logFC) >=  log2(2))%>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)%>%filter(zfin_id_symbol %in% epi_reg)
  name=paste0("epi_reg_cpm.",i)
  assign(name,deg)
}

for (i in c(5:8)) {
  deg = get(paste0("LD.all.",i))
  deg = deg %>% filter(mean.cpm.pos > 5|mean.cpm.wt >5)%>% filter(t.test < 0.05)%>%filter(FDR < 0.01)%>%
    filter(abs(logFC) >=  log2(2))%>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)%>%filter(zfin_id_symbol %in% DNA_modi)
  name=paste0("DNA_modi_cpm.",i)
  assign(name,deg)
}

for (i in c(5:8)) {
  deg = get(paste0("LD.all.",i))
  deg = deg %>% filter(mean.cpm.pos > 5|mean.cpm.wt >5)%>% filter(t.test < 0.05)%>%filter(FDR < 0.01)%>%
    filter(abs(logFC) >=  log2(2))%>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)%>%filter(zfin_id_symbol %in% hist_modi)
  name=paste0("hist_modi_cpm.",i)
  assign(name,deg)
}

for (i in c(5:8)) {
  deg = get(paste0("LD.all.",i))
  deg = deg %>% filter(mean.cpm.pos > 5|mean.cpm.wt >5)%>% filter(t.test < 0.05)%>%filter(FDR < 0.01)%>%
    filter(abs(logFC) >=  log2(2))%>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)%>%filter(zfin_id_symbol %in% RNA_modi)
  name=paste0("RNA_modi_cpm.",i)
  assign(name,deg)
}


for (i in c(5:8)) {
  deg = get(paste0("LD.all.",i))
  deg = deg %>% filter(mean.cpm.pos > 5|mean.cpm.wt >5)%>% filter(t.test < 0.05)%>%filter(FDR < 0.01)%>%
    filter(abs(logFC) >=  log2(2))%>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)%>%filter(zfin_id_symbol %in% RNA_proc)
  name=paste0("RNA_proc_cpm.",i)
  assign(name,deg)
} 

#trans_temporal GC_affected_sig. epi
m.epi_reg = intersect(unique(c(epi_reg_cpm.5$zfin_id_symbol,epi_reg_cpm.6$zfin_id_symbol)),
                       unique(c(epi_reg_cpm.7$zfin_id_symbol,epi_reg_cpm.8$zfin_id_symbol)))

m.dna.modi = intersect(unique(c(DNA_modi_cpm.5$zfin_id_symbol,DNA_modi_cpm.6$zfin_id_symbol)),
                         unique(c(DNA_modi_cpm.7$zfin_id_symbol,DNA_modi_cpm.8$zfin_id_symbol)))

m.hist_mod = intersect(unique(c(hist_modi_cpm.5$zfin_id_symbol,hist_modi_cpm.6$zfin_id_symbol)),
                       unique(c(hist_modi_cpm.7$zfin_id_symbol,hist_modi_cpm.8$zfin_id_symbol)))

m.RNA.modi = intersect(unique(c(RNA_modi_cpm.5$zfin_id_symbol,RNA_modi_cpm.6$zfin_id_symbol)),
                       unique(c(RNA_modi_cpm.7$zfin_id_symbol,RNA_modi_cpm.8$zfin_id_symbol)))

m.RNA_proc = intersect(unique(c(RNA_proc_cpm.5$zfin_id_symbol,RNA_proc_cpm.6$zfin_id_symbol)),
                       unique(c(RNA_proc_cpm.7$zfin_id_symbol,RNA_proc_cpm.8$zfin_id_symbol)))
#all_sig_early

all_sig_early_epi = rbind(DNA_modi_cpm.5,
                         hist_modi_cpm.5,hist_modi_cpm.6,
                         RNA_modi_cpm.5,RNA_modi_cpm.6,
                         RNA_proc_cpm.5,RNA_proc_cpm.6)


all_sig_early_epi[,"epi"]=c(paste0(DNA_modi_cpm.5$zfin_id_symbol," (DNA_modi.)"),
                            paste0(hist_modi_cpm.5$zfin_id_symbol," (hist_modi.)"),
                            paste0(hist_modi_cpm.6$zfin_id_symbol," (hist_modi.)"),
                            paste0(RNA_modi_cpm.5$zfin_id_symbol," (RNA_modi.)"),
                            paste0(RNA_modi_cpm.6$zfin_id_symbol," (RNA_modi.)"),
                            paste0(RNA_proc_cpm.5$zfin_id_symbol," (RNA_proc.)"),
                            paste0(RNA_proc_cpm.6$zfin_id_symbol," (RNA_proc.)"))
#all_sig_prepostLD

all_sig_ppLD_epi = rbind(DNA_modi_cpm.7,DNA_modi_cpm.8,
                         hist_modi_cpm.7,hist_modi_cpm.8,
                         RNA_modi_cpm.7,RNA_modi_cpm.8,
                         RNA_proc_cpm.7,RNA_proc_cpm.8)

all_sig_ppLD_epi[,"epi"]=c( paste0(DNA_modi_cpm.7$zfin_id_symbol," (DNA_modi.)"),
                            paste0(DNA_modi_cpm.8$zfin_id_symbol," (DNA_modi.)"),
                            paste0(hist_modi_cpm.7$zfin_id_symbol," (hist_modi.)"),
                            paste0(hist_modi_cpm.8$zfin_id_symbol," (hist_modi.)"),
                            paste0(RNA_modi_cpm.7$zfin_id_symbol," (RNA_modi.)"),
                           paste0(RNA_modi_cpm.8$zfin_id_symbol," (RNA_modi.)"),
                           paste0(RNA_proc_cpm.7$zfin_id_symbol," (RNA_proc.)"),
                            paste0(RNA_proc_cpm.8$zfin_id_symbol," (RNA_proc.)"))
LDgene_epi = rbind(all_sig_early_epi,all_sig_ppLD_epi)
write.csv(LDgene_epi,"./outputs/s.table.9_LDgene_epi.csv")
LDgene_epi_dict =strsplit(unique(LDgene_epi[,"epi"]),split = " ")

LDgene_epi= c()
for (i in 1:75) {
  gn =  LDgene_epi_dict[[i]][1]
  LDgene_epi=c(LDgene_epi,gn)
}
LDgene_epi_type= c()
for (i in 1:75) {
  type =  LDgene_epi_dict[[i]][2]
  LDgene_epi_type=c(LDgene_epi_type,type)
}

names(LDgene_epi_type)=LDgene_epi
#write.csv(LDgene_epi,"./outputs/s.table.10_LDgene_epi.csv")

all_sig_epi = unique(c(all_sig_early_epi$zfin_id_symbol,all_sig_ppLD_epi$zfin_id_symbol))
all_sig_epi_cpm = filter(LD.all.8,zfin_id_symbol %in% all_sig_epi)


#map for epi

all_sig_epi_cpm_table = data.frame()

for(i in c(5:8)){
  deg = get(paste0("deg.cpm",i))
  deg = filter(deg, zfin_id_symbol %in% all_sig_epi)
  all_sig_epi_cpm_table = rbind(deg, all_sig_epi_cpm_table)
}

#matrix for heatmap

sig_epi_ht = as.matrix(data.frame( "d6" = all_sig_epi_cpm_table$logFC[c(which(all_sig_epi_cpm_table$source=="d6"))],
                       "d13"= all_sig_epi_cpm_table$logFC[c(which(all_sig_epi_cpm_table$source=="d13"))],
                       "d120"= all_sig_epi_cpm_table$logFC[c(which(all_sig_epi_cpm_table$source=="d120"))],
                       "LD"= all_sig_epi_cpm_table$logFC[c(which(all_sig_epi_cpm_table$source=="LD"))]))
sig_epi_ht.FDR = as.matrix(data.frame( "d6" = all_sig_epi_cpm_table$FDR[c(which(all_sig_epi_cpm_table$source=="d6"))],
                                   "d13"= all_sig_epi_cpm_table$FDR[c(which(all_sig_epi_cpm_table$source=="d13"))],
                                   "d120"= all_sig_epi_cpm_table$FDR[c(which(all_sig_epi_cpm_table$source=="d120"))],
                                   "LD"= all_sig_epi_cpm_table$FDR[c(which(all_sig_epi_cpm_table$source=="LD"))]))


sig_epi_ht.FDR[sig_epi_ht.FDR <= 0.01] <-"*"
sig_epi_ht.FDR[sig_epi_ht.FDR > 0.01] <-NA
sig_epi_ht.FDR[abs(sig_epi_ht) < log2(1.5)] <-NA



st_sig_tx_1= sig_epi_ht
distfunc <- function(x) dist(x, method="euclidean")
hclustfunc <- function(x) hclust(x, method="mcquitty")

r.ns = all_sig_epi_cpm_table$zfin_id_symbol[c(which(all_sig_epi_cpm_table$source=="d6"))]
row.names(st_sig_tx_1)= paste(r.ns,LDgene_epi_type[r.ns])

hc = hclustfunc(distfunc(st_sig_tx_1))
ht.label.order = hc$labels[hc$order]
write(ht.label.order,"./figures/Fig6.epi_LD.label.order.txt")

rhc = hclust(dist(st_sig_tx_1, method="euclidean"), method="mcquitty")
gr.row <- cutree(rhc, 5)
tree.table = data.frame("gene" = names(gr.row), "cluster"= gr.row)
write.csv(tree.table,"./figures/Fig6.epi_LD..cluster.order.csv")
library(viridis)
desturate = viridis(3)
row.col = LDgene_epi_type[r.ns]
row.col[row.col == "(DNA_modi.)"] <-1
row.col[row.col == "(hist_modi.)"] <-2
row.col[row.col == "(RNA_modi.)"] <-3
row.col[row.col == "(RNA_proc.)"] <-3


tiff(paste0("./figures/Fig6a.epi.LD.DEGS.2.heatmap.",nrow(st_sig_tx_1) ,".tiff"), height = 15,width = 8, res=1200, units = "cm", compression = "lzw")

try(heatmap.2((st_sig_tx_1), col=mycol, scale='none',
              dendrogram = "row",Colv = F,breaks=breaks,hclustfun=hclustfunc, distfun=distfunc,
              density='density', trace = 'none',
              na.color = F,
              offsetRow = 0,
              offsetCol = 0,
              cexCol = 1,
              cexRow = 1,
              #colsep = 4,
              labCol = c("6","13","120","LD"),
              RowSideColors=desturate[as.numeric(row.col)],
              #lwid = c(1,2.5),
              cellnote=sig_epi_ht.FDR,
              notecex=0.5,
              notecol="black",
              labRow = "",
              margins = c(5, 15),
              keysize = 0.5 ), silent = T)

dev.off()


tiff(paste0("./figures/Fig6a.epi.LD.DEGS.2.labeled.heatmap.",nrow(st_sig_tx_1) ,".tiff"), height = 30,width = 10, res=1200, units = "cm", compression = "lzw")

try(heatmap.2((st_sig_tx_1), col=mycol, scale='none',
              dendrogram = "row",Colv = F,breaks=breaks,hclustfun=hclustfunc, distfun=distfunc,
              density='density', trace = 'none',
              na.color = F,
              offsetRow = 0,
              offsetCol = 0,
              cexCol = 1,
              cexRow = 1,
              #colsep = 4,
              labCol = c("6","13","120","LD"),
              RowSideColors=desturate[as.numeric(row.col)],
              #lwid = c(1,2.5),
              cellnote=sig_epi_ht.FDR,
              notecex=0.5,
              notecol="black",
              #labRow = "",
              margins = c(5, 15),
              keysize = 0.5 ), silent = T)

dev.off()

epi.bar.table = data.frame("gene" = LDgene_epi,"type" = LDgene_epi_type)

epi.bar= count(epi.bar.table, type)[c(1:3),]
epi.bar[1,1] = "DNA" #3
epi.bar[2,1] = "Histone" #22
epi.bar[3,1] = "RNA"
epi.bar[3,2]=50 # 6(from RNA modification) + 44(from RNA_process) - 5 overlapped RNA modi and proc.

tiff(paste0("./figures/Fig6.epi.LD.DEGS.2.bar.",nrow(epi.bar.table) ,".tiff"), height = 7,width = 6, res=1200, units = "cm", compression = "lzw")

epi.bar[,"modifications"]= "modifications"
ggplot(epi.bar, aes(x=modifications, y= n, fill=type)) +
  geom_bar(position="stack",stat="identity") +
  scale_fill_viridis(discrete=TRUE, name="") +
  theme_classic() +
  ylab("") +
  xlab("")
dev.off()



###s.fig7 # comparison with PNAS, https://www.pnas.org/doi/10.1073/pnas.1820842116, Provençal et al., 2020
#Supplementary table 5: List of significantly regulated Pro-diff+WO +acute transcripts that map to a long-lasing DMS (n=702).
s5_pnas = read.csv("./supplementary/pnas.1820842116.sd01_st5.csv")


#LD
ft.LD.deg = deg.cpm8 %>% filter(t.test <0.05)%>% filter(FDR <0.05) %>% filter(mean.cpm.pos > 5|mean.cpm.wt >5)

hz_LD_list.up = filter(zh.all,zfin_id_symbol %in% ft.LD.deg$zfin_id_symbol[c(which(ft.LD.deg$logFC > log2(1.3)))])
hz_LD_list.dn = filter(zh.all,zfin_id_symbol %in% ft.LD.deg$zfin_id_symbol[c(which(ft.LD.deg$logFC < -log2(1.3)))])


s5_intLD.up= s5_pnas %>% filter(FC_pro.diff.WO.acute >0) %>% 
  filter(Gene_symbol%in% unique(hz_LD_list.up$human.Symbol))
s5_intLD.dn= s5_pnas %>% filter(FC_pro.diff.WO.acute <0) %>% 
  filter(Gene_symbol%in% unique(hz_LD_list.dn$human.Symbol))

#13
ft.13.deg = deg.cpm6 %>% filter(t.test <0.05)%>% filter(FDR <0.05) %>%   filter(mean.cpm.pos > 5|mean.cpm.wt >5)

hz_13_list.up = filter(zh.all,zfin_id_symbol %in% 
                                     ft.13.deg$zfin_id_symbol[c(which(ft.13.deg$logFC > log2(1.3)))])
hz_13_list.dn = filter(zh.all,zfin_id_symbol %in% 
                                     ft.13.deg$zfin_id_symbol[c(which(ft.13.deg$logFC < -log2(1.3)))])


s5_intd13.up= s5_pnas %>% filter(FC_pro.diff.WO.acute >0) %>% 
  filter(Gene_symbol%in% hz_13_list.up$human.Symbol)
s5_intd13.dn= s5_pnas %>% filter(FC_pro.diff.WO.acute <0) %>% 
  filter(Gene_symbol%in% hz_13_list.dn$human.Symbol)


s5_up_list = list(
  "overlap_d13" = unique(s5_intd13.up$Gene_symbol),
  "overlap_LD" = unique(s5_intLD.up$Gene_symbol)
)

s5_dn_list = list(
  "overlap_d13" = unique(s5_intd13.dn$Gene_symbol),
  "overlap_LD" = unique(s5_intLD.dn$Gene_symbol)
)



library(eulerr)

    
fit.up= euler(s5_up_list)
fit.dn= euler(s5_dn_list, shape="ellipse")


venncol6 = c("#A50026","#F46D43","#FEE090" ,"#313695","#4575B4","#ABD9E9")

up.venn =plot(fit.up, 
              quantities = TRUE,
              fill = venncol6[1:2], alpha =0.3,
              lty = 1,
              labels = list(font = 4))

dn.venn =plot(fit.dn, 
              quantities = TRUE,
              fill = venncol6[6:5], alpha =0.3,
              lty = 1,
              labels = list(font = 4))


tiff("./figures/s.Fig7_deg.all.venn3.15x12.tiff", width =20 , height = 10, res = 1200, units = "cm",compression = "lzw")
#print(deg.all.venn)
plot_grid("",up.venn,"",dn.venn,"",nrow = 1,label_size = 12,rel_widths =c(0.2,1,0.4,1,0.2) )
dev.off()









