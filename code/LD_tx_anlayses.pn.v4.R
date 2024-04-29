#####s.pos vs s.neg
#selection of DEG established in early life


#DEG after LD between star:bPAC+ vs bPAC- 
sig7 = deg.cpm12 %>% filter(FDR < 0.05) %>% filter(mean.cpm.sub > 5|mean.cpm.cont >5)%>% 
  filter(abs(logFC) >= log2(1.5))%>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)

#DEG after LD in star:bPAC+
sig8 = filter(deg.cpm13, FDR < 0.05)%>% filter(mean.cpm.sub > 5|mean.cpm.cont >5)%>% 
  filter(abs(logFC) >=  log2(1.5))%>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)

#non-deg after LD in star:bPAC+
sig8.non = filter(deg.cpm13, FDR > 0.05) %>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)%>% filter(mean.cpm.sub > 5|mean.cpm.cont >5)

#DEG after LD in bPAC-
sig9 = filter(deg.cpm15, FDR < 0.05)%>% filter(mean.cpm.sub > 5|mean.cpm.cont >5)%>% 
  filter(abs(logFC) >=  log2(1.5))%>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)

#Common DEGs following LD
#DEG pre LD in bPAC-
sig10 = filter(deg.cpm11, FDR < 0.05)%>% filter(mean.cpm.sub > 5|mean.cpm.cont >5)%>% 
  filter(abs(logFC) >=  log2(1.5))%>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)

#overlapping bPAC- vs.bpac-LD between bPAC vs. bPAC+LD
sig89_genes=intersect(sig8$zfin_id_symbol,sig9$zfin_id_symbol)#85 genes

#check direction of FC
sig8$zfin_id_symbol[c(which(sig8[c(which(sig8$zfin_id_symbol%in%sig89_genes)),"logFC"]*sig9[c(which(sig9$zfin_id_symbol%in%sig89_genes)),"logFC"] <0))]
#if it's 0. proceed to next

sig89_8 = sig8[c(which(sig8$zfin_id_symbol%in%sig89_genes)),]
sig89_9 = sig9[c(which(sig9$zfin_id_symbol%in%sig89_genes)),]
sig89 = rbind(sig89_8[,-c(subj.cpm.list.LD[[1]],cont.cpm.list.LD[[1]])],sig89_9[,-c(subj.cpm.list.LD[[3]],cont.cpm.list.LD[[3]])])

#check overlapping with LD
deg.cpm12[c(which(deg.cpm12$zfin_id_symbol%in%sig89_genes)),] %>% 
  filter(abs(logFC)>log2(1.5)) %>% filter(FDR <0.05)
#if it's 0. proceed to next, 

#.. overlapping 17 genes show same direction of regulation, and are not part of post-LD-degs 
#.. This result indicates that those DEGs are involved in normal response to LD.

# DEG2 = bPAC+ LD vs. bPAC+ d120, DEG4 = bPAC- LD vs bPAC- d120, DEG5 = bPAC+ LD vs. bPAC- LD
deg.list2 = list(
  "DEG 2" = sig8$zfin_id_symbol,
  "DEG 4" = sig9$zfin_id_symbol
  )

#check direction of FC
sig79_genes=intersect(sig7$zfin_id_symbol,sig9$zfin_id_symbol) # 34 genes
sig79_genes_only=sig79_genes[-c(which(sig79_genes %in% intersect(sig79_genes,sig8$zfin_id_symbol)))]
sig9$zfin_id_symbol[c(which(sig9[c(which(sig9$zfin_id_symbol%in%sig79_genes)),"logFC"]*sig9[c(which(sig7$zfin_id_symbol%in%sig79_genes)),"logFC"] <0))]
#.. DEG4 and DEG6 have 34 overlapping gene but oppositely regulated except egr1. Thus, they are part of post-LD_DEGs except egr1


deg24= unique(c(deg.list2$`DEG 4`,deg.list2$`DEG 2`))
deg24.primed = deg24[-c(which(deg24 %in% c(intersect(deg.list2$`DEG 4`,deg.list2$`DEG 2`))))]
deg24.primed = deg24.primed[-c(which(deg24.primed %in% sig79_genes_only))]
write.csv(deg24.primed, "./outputs/s.table3_list of deg24.primed.pn.csv")

deg245= unique(c(deg24,deg.list2$`DEG 5`))
deg245.primed= deg245[-c(which(deg24 %in% c(intersect(deg.list2$`DEG 4`,deg.list2$`DEG 2`))))]
write.csv(deg245.primed, "./outputs/s.table3_list of primed_plus_DEG5.pn.csv")



ab.LD.res = unique(c(sig8$id,sig7$id,sig9$id))
ab.LD.res=ab.LD.res[- c(which(ab.LD.res %in% unique(sig89$id)))]
ab.LD.res.deg.gene = intersect(deg.cpm9[c(which(deg.cpm9$id %in% ab.LD.res )),8],prot_id_map$SYMBOL)

LD.deg.list = rbind(sig7[,-c(subj.cpm.list[[12]],cont.cpm.list[[12]])],
                    sig8[,-c(subj.cpm.list.LD[[1]],cont.cpm.list.LD[[1]])],
                    sig9[,-c(subj.cpm.list.LD[[3]],cont.cpm.list.LD[[3]])])
LD.deg.list = LD.deg.list[-c(which(LD.deg.list$zfin_id_symbol %in% unique(sig89$zfin_id_symbol))),]
LD.deg.list$source= gsub("sp.m4-LD",x = LD.deg.list$source ,"DEG2")
LD.deg.list$source= gsub("sn.m4-LD",x = LD.deg.list$source ,"DEG1")
LD.deg.list$source= gsub("LDpn",x = LD.deg.list$source ,"DEG3")

LD.deg.list
#s.table3_list of LD-DEGs
LD.deg.list.table=LD.deg.list
write.csv(LD.deg.list.table, "./outputs/s.table3_list of LD-DEGs.pn.csv")


###venn diagram
'
deg.list2 = list(
  "DEG 2" = sig9$zfin_id_symbol,
  "DEG 4" = sig8$zfin_id_symbol
  )
'

desturate=viridis(3)[c(2,1)]

library(eulerr)

set.seed(7)#99
fit3= euler(deg.list2,shape = "ellipse")
fit3
deg.venn =  plot(fit3, 
                 quantities = T,
                 fill = desturate, alpha =0.3,
                 lty = 1,
                 labels = ""
)

#

'                   original fitted residuals regionError
DEG 2                  545    545         0           0
DEG 4                   14     14         0           0
DEG 5                  488    488         0           0
DEG 2&DEG 4             54     54         0           0
DEG 2&DEG 5           1124   1124         0           0
DEG 4&DEG 5              3      3         0           0
DEG 2&DEG 4&DEG 5       31     31         0           0

diagError: 0 
stress:    0

'


tiff("./figures/Fig5E_venn.LD.deg.12x12.pn_wo_labela.tiff",  height = 8,width = 8, res=1200, units = "cm", compression = "lzw")
plot(deg.venn)
dev.off()



#####
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

LD.deg_all.up=LD.deg.list%>%filter(logFC > 0)%>%filter(zfin_id_symbol %in% deg24.primed)
zh.LD_DEG.up= unique(filter(zh.all, zfin_id_symbol %in% LD.deg_all.up$zfin_id_symbol)[,2])


#up-regulated
LD.deg_all.up=LD.deg.list%>%filter(logFC > 0)%>%filter(zfin_id_symbol %in% deg24.primed)
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
write.csv(chea3_LD.int.mean.rank,"./outputs/s.table5_chea3_LD_DEG.up_results_Integrated--meanRank.pn.csv")  


chea3 = data.frame(chea3_LD$`Integrated--meanRank`)
#chea3.sb50= chea3 %>% filter(as.numeric(Rank) <= 50) 
chea3.up.deg = chea3 %>% filter(as.numeric(Rank) <= 50) %>%
  filter(TF %in% c(zh.all$human.Symbol[c(which(zh.all$zfin_id_symbol %in% deg_all.up$zfin_id_symbol))]))

library(forcats) 

top.tf.LD=chea3.up.deg$TF
top.tf.LD.rank= c()
for (i in 1:length(chea3.dn.deg$Query.Name)) {
  counts=c(top.tf.LD.rank,length(str_split(chea3.dn.deg$Overlapping_Genes,pattern = ",")[[i]]))
  top.tf.LD.rank = c(top.tf.LD.rank,counts)
}
top.tf.LD.rank= paste0("(",top.tf.LD.rank,")")
top.tf.LD.score= log2(1632/as.numeric(chea3.up.deg$Score))

#.. is part of LD_DEG? 
topTF.zh_LD=filter(LD.deg.list, zfin_id_symbol %in% unique(c(filter(zh.all,human.Symbol %in% top.tf.LD)$zfin_id_symbol)))$zfin_id_symbol
topTF.hz_LD=unique(c(filter(zh.all,zfin_id_symbol %in% topTF.zh_LD)$human.Symbol))

#.. "CSRNP1" "EGR1"   "FOSL2"  "JUNB"   "EN2"    "FEZF2"  "FOS"   ,for up-reg LD-dges.  
#.. "FEZF2"  "CSRNP1" "FOS"    "BATF2"  are part of LD_DEGs

#top10 TF .but only 9 are detected
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


tiff("./figures/s.Fig6a_LD-up.tf.pn.tiff", width = 12, height = 12,units = "cm", res = 1200, compression = "lzw")
plot(LD.tf.p)
dev.off()

#down-regulated
LD.deg_all.dn=LD.deg.list%>%filter(logFC < 0)%>%filter(zfin_id_symbol %in% deg24.primed)

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
write.csv(chea3_LD.int.mean.rank,"./outputs/s.table5_chea3_LD_DEG.dn_results_Integrated--meanRank.pn.csv")  

chea3 = data.frame(chea3_LD$`Integrated--meanRank`)
#chea3.sb50= chea3 %>% filter(as.numeric(Rank) <= 50) 
chea3.dn.deg = chea3 %>% filter(as.numeric(Rank) <= 50) %>%
  filter(TF %in% c(zh.all$human.Symbol[c(which(zh.all$zfin_id_symbol %in% deg_all.dn$zfin_id_symbol))]))

library(forcats) 

top.tf.LD=chea3.dn.deg$TF
top.tf.LD.rank= c()
for (i in 1:length(chea3.dn.deg$Query.Name)) {
  counts=c(top.tf.LD.rank,length(str_split(chea3.dn.deg$Overlapping_Genes,pattern = ",")[[i]]))
  top.tf.LD.rank = c(top.tf.LD.rank,counts)
  }
top.tf.LD.rank= paste0("(",top.tf.LD.rank,")")
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


tiff("./figures/s.Fig6b_LD-down.tf.pn.tiff", width = 12, height = 12,units = "cm", res = 1200, compression = "lzw")
plot(LD.tf.p)
dev.off()



# GO enrichment analysis of LD-DEGs ---------------------------------------


#list of background genes
for (i in c(12,13:15)){
  deg = get(paste0("deg.cpm",i))
  deg = deg %>% filter(mean.cpm.sub > 5|mean.cpm.cont >5)
  deg = prot_id_map[deg$zfin_id_symbol, "ENTREZID"]
  deg = na.omit(deg)
  name = paste0("bk.gene.LD",i)
  assign(name, deg)
  print(i)
}

bk.all.LD = unique(c(bk.gene.LD12,bk.gene.LD13,bk.gene.LD15))

#each time point
time.name2=c("DEG","DEG2","DEG4","DEG5")

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

write.csv(LD.tc_enrich.up.tb,"./outputs/s.table4a_LD.tc_enrich.up.pn.tb.csv" )
write.csv(LD.tc_enrich.dn.tb,"./outputs/s.table4b_LD.tc_enrich.dn.pn.tb.csv" )



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
tryCatch({
  for (k in 1:3) {
    
    ##up
    tc_enrich.up <- clusterProfiler::compareCluster(tc_point_up, fun = "enrichGO",  OrgDb='org.Dr.eg.db', ont = gos[k],
                                                    bk.all.LD,pvalueCutoff = 0.05,
                                                    pAdjustMethod = "BH",qvalueCutoff = 0.2,
                                                    minGSSize = 10,
                                                    maxGSSize = 600, 
                                                    readable = TRUE)
    
    tc_enrich.up.tb = data.frame(tc_enrich.up[1:length(tc_enrich.up[]$ID)])
    tc_enrich.up.tb[,"source"]<-gos[k]
    all.LD.tc_enrich.up.tb =rbind(all.LD.tc_enrich.up.tb,tc_enrich.up.tb)
    
    
  }
}, error=function(e){})


for (k in 1:3) {  
  ##down
  tc_enrich.dn <- clusterProfiler::compareCluster(tc_point_dn, fun = "enrichGO",  OrgDb='org.Dr.eg.db', ont = gos[k],
                                                  bk.all.LD,pvalueCutoff = 0.05,
                                                  pAdjustMethod = "BH",qvalueCutoff = 0.2,
                                                  minGSSize = 10,
                                                  maxGSSize = 1000, 
                                                  readable = TRUE)
  
  tc_enrich.dn.tb = data.frame(tc_enrich.dn[1:length(tc_enrich.dn[]$ID)])
  tc_enrich.dn.tb[,"source"]=gos[k]
  
  
  
  all.LD.tc_enrich.dn.tb =rbind(all.LD.tc_enrich.dn.tb, tc_enrich.dn.tb) 
} 


tryCatch({
  
  a= str_split(all.LD.tc_enrich.up.tb[,4],pattern = "/")
  if (length(a)>1) {
    b= data.frame(Reduce(rbind,a))
  } else {b= data.frame(a[[1]][1],a[[1]][2])  
  }  
  all.LD.tc_enrich.up.tb[,"n.deg"] = as.numeric(b[,2])
  
  a= str_split(all.LD.tc_enrich.up.tb[,5],pattern = "/")
  if (length(a)>1) {
    b= data.frame(Reduce(rbind,a))
  } else {b= data.frame(a[[1]][1],a[[1]][2])  
  }  
  all.LD.tc_enrich.up.tb[,"gene.in.catg"] = as.numeric(b[,1])
  
  all.LD.tc_enrich.up.tb[,"ratio.catg"] = all.LD.tc_enrich.up.tb[,10]/all.LD.tc_enrich.up.tb[,"gene.in.catg"]
  
  #filtering
  #all.LD.tc_enrich.up.tb= all.LD.tc_enrich.up.tb %>% filter(Count >= 10) %>% filter(ratio.catg >= 0.05)
},error=function(e){})


a= str_split(all.LD.tc_enrich.dn.tb[,4],pattern = "/")
if (length(a)>1) {
  b= data.frame(Reduce(rbind,a))
} else {b= data.frame(a[[1]][1],a[[1]][2])  
}
all.LD.tc_enrich.dn.tb[,"n.deg"] = as.numeric(b[,2])

a= str_split(all.LD.tc_enrich.dn.tb[,5],pattern = "/")
if (length(a)>1) {
  b= data.frame(Reduce(rbind,a))
} else {b= data.frame(a[[1]][1],a[[1]][2])  
}
all.LD.tc_enrich.dn.tb[,"gene.in.catg"] = as.numeric(b[,1])

all.LD.tc_enrich.dn.tb[,"ratio.catg"] = c(all.LD.tc_enrich.dn.tb[,10]/all.LD.tc_enrich.dn.tb[,"gene.in.catg"])

#filtering
all.LD.tc_enrich.dn.tb= all.LD.tc_enrich.dn.tb %>% filter(Count >= 10) %>% filter(ratio.catg >= 0.05)



#write_combined_table

write.csv(all.LD.tc_enrich.dn.tb,"./outputs/s.table4d_all_LD.tc_enrich.dn.pn.tb.csv" )


##### GO-primed

#read table 
prime24_GO = read.csv("./data/all_primed genes_gProfiler_drerio_02-08-2023_15-33-14__intersections.csv")

prime24_GO.genes = list(
  #top CC
  "synapse"=str_split(prime24_GO[c(which(prime24_GO$term_name == "synapse")),11],pattern = ",")[[1]],
  #top BP
  "cell_morpho"=str_split(prime24_GO[c(which(prime24_GO$term_name == "cell morphogenesis")),11],pattern = ",")[[1]],
  "reg_signaling"=str_split(prime24_GO[c(which(prime24_GO$term_name == "regulation of signaling")),11],pattern = ",")[[1]],
  #top MF
  "salt.tta"=str_split(prime24_GO[c(which(prime24_GO$term_name == "salt transmembrane transporter activity")),11],pattern = ",")[[1]],
  "inorganic.itt"=str_split(prime24_GO[c(which(prime24_GO$term_name == "inorganic ion transmembrane transport")),11],pattern = ",")[[1]]
  
)


#primed.synapse+cell_morph
synapse.cell_morpho=intersect(prime24_GO.genes[[1]],prime24_GO.genes[[2]])
#primed.synapse+reg_signaling
synapse.reg_sig=intersect(prime24_GO.genes[[1]],prime24_GO.genes[[3]])
#primed.synapse+salt.transmembrane transporter activity
synapse.salt.tta=intersect(prime24_GO.genes[[1]],prime24_GO.genes[[4]])


library(stringr)
deg.cpm12[c(which(deg.cpm12$zfin_id_symbol %in% str_to_lower(synapse.cell_morpho))),c(1,5,8,18,19)]
deg.cpm12[c(which(deg.cpm12$zfin_id_symbol %in% str_to_lower(synapse.reg_sig))),c(1,5,8,18,19)]






######Disgenet
# disease associated_enrichment_disignet ----------------------------------
library(disgenet2r)
#disgenet_api_key <- get_disgenet_api_key(
#  email = "user@email.com", 
#  password = "myspwd" )
#Sys.setenv(DISGENET_API_KEY= disgenet_api_key)

#

# for GC-altered
zh_LD_list.deg6= unique(filter(zh.all, zfin_id_symbol %in% c(deg.list2$`DEG 5`))[,2])

try(res_enrich <-disease_enrichment( entities =zh_LD_list.deg6, vocabulary = "HGNC",
                                     database = "CURATED" ) # "deg5"
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

write.csv(table.dis.LD.sig,"./outputs/s.table7_DisGeNET_enrichment_deg6.csv")

tiff("./figures/Fig6A_top.disgenet.deg6.tiff", width = 13, height = 12,units = "cm", res = 1200, compression = "lzw")
plot(res_enrich, class = "Enrichment", count =20,  cutoff= 0.05, nchars=70)
dev.off()


# for GC-primed genes
zh_LD_list= unique(filter(zh.all, zfin_id_symbol %in% deg24.primed)[,2])

try(res_enrich <-disease_enrichment( entities =zh_LD_list, vocabulary = "HGNC",
                                     database = "CURATED" ) # "DEG2+4"
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

write.csv(table.dis.LD.sig,"./outputs/s.table6_DisGeNET_enrichment_deg24.csv")

tiff("./figures/Fig6A_top.disgenet_deg24.tiff", width = 13, height = 10,units = "cm", res = 1200, compression = "lzw")
plot(res_enrich, class = "Enrichment", count =40,  cutoff= 0.05, nchars=70)
dev.off()


#Venn diagram
dis.name= unique(table.dis.LD.sig$Description)

schiz.genes = unique(c(unlist(str_split(table.dis.LD.sig$shared_symbol[c(which(table.dis.LD.sig$Description %in% dis.name[c(1)]))],pattern = ";"))) ) 
bipolar.genes = unique(c(unlist(str_split(table.dis.LD.sig$shared_symbol[c(which(table.dis.LD.sig$Description %in% dis.name[c(2)]))],pattern = ";"))) ) 
int.dis.genes = unique(c(unlist(str_split(table.dis.LD.sig$shared_symbol[c(which(table.dis.LD.sig$Description %in% dis.name[c(3)]))],pattern = ";"))) ) 
auti.genes = unique(c(unlist(str_split(table.dis.LD.sig$shared_symbol[c(which(table.dis.LD.sig$Description %in% dis.name[c(4)]))],pattern = ";"))) ) 
dep.genes = unique(c(unlist(str_split(table.dis.LD.sig$shared_symbol[c(which(table.dis.LD.sig$Description %in% dis.name[c(5,7)]))],pattern = ";"))) )

schiz.disig.zh = zh.all %>% filter(human.Symbol %in% schiz.genes)
bipolar.disig.zh = zh.all %>% filter(human.Symbol %in% bipolar.genes)
int.dis.disig.zh = zh.all %>% filter(human.Symbol %in% int.dis.genes)
auti.disig.zh = zh.all %>% filter(human.Symbol %in% auti.genes)
dep.disig.zh = zh.all %>% filter(human.Symbol %in% dep.genes)

schiz.disig.LD = deg.cpm12 %>% filter(zfin_id_symbol %in% intersect(ab.LD.res.deg.gene,unique(schiz.disig.zh$zfin_id_symbol)))
bipolar.disig.LD = deg.cpm12 %>% filter(zfin_id_symbol %in% intersect(ab.LD.res.deg.gene,unique(bipolar.disig.zh$zfin_id_symbol)))
int.dis.disig.LD = deg.cpm12 %>% filter(zfin_id_symbol %in% intersect(ab.LD.res.deg.gene,unique(int.dis.disig.zh$zfin_id_symbol)))
auti.disig.LD = deg.cpm12 %>% filter(zfin_id_symbol %in% intersect(ab.LD.res.deg.gene,unique(auti.disig.zh$zfin_id_symbol)))
dep.disig.LD = deg.cpm12 %>% filter(zfin_id_symbol %in% intersect(ab.LD.res.deg.gene,unique(dep.disig.zh$zfin_id_symbol))) 

###

disig.list = list(
  "schizo." = unique(schiz.disig.LD$zfin_id_symbol),
  "bipolar" = unique(bipolar.disig.LD$zfin_id_symbol),
  "int.dis" = unique(int.dis.disig.LD$zfin_id_symbol),
  "Autistic" = unique(auti.disig.LD$zfin_id_symbol),
  "depressive" = unique(dep.disig.LD$zfin_id_symbol)
  
)


library(eulerr)


fit.dis= euler(disig.list) # 10 genes between dep. and schizo. are missing because of the adjustment. manually adjust.
venncol4 = viridis(5)

disnet.all.venn =plot(fit.dis, 
                      quantities = T,
                      fill = venncol4, alpha =0.3,
                      lty = 1,
                      labels = "" #list(font = 4)
)


# need to add losses in fitted  depressive&schizo-> 10
tiff("./figures/Fig6B_disgnet.primed.venn5.pn.12x12.nolabel.tiff", height = 12,width = 12, res=1200, units = "cm", compression = "lzw")
plot(disnet.all.venn)
dev.off()


##### version for ALL DEGs

#selection
for (i in c(9:12)) {
  deg = get(paste0("deg.cpm",i))
  deg = deg %>% filter(zfin_id_symbol %in% deg245.primed)
  name=paste0("LD.all.",i)
  assign(name,deg)
}
deg245.primed

all_LD_deg=rbind(LD.all.9[,-c(subj.cpm.list[[9]],cont.cpm.list[[9]])],
                 LD.all.10[,-c(subj.cpm.list[[10]],cont.cpm.list[[10]])],
                 LD.all.11[,-c(subj.cpm.list[[11]],cont.cpm.list[[11]])],
                 LD.all.12[,-c(subj.cpm.list[[12]],cont.cpm.list[[12]])])
write.csv(all_dis_LD_deg,"./outputs/tt_asp_LD_deg.pn.csv")

#trans_temporal, DEG in LD and d6 + (d120 or d13)
main.LD.all.FDR = Reduce(intersect, list(LD.all.9$zfin_id_symbol[c(which(LD.all.9$FDR <0.05))],
                                         LD.all.12$zfin_id_symbol[c(which(LD.all.12$FDR <0.05))]))

matched.alt = c()

for (i in 1:nrow(LD.all.12)) {
  if (abs(LD.all.12$logFC[i])> log2(1.5) &
       abs(LD.all.9$logFC[i])> log2(1.5)){
    if (LD.all.12$logFC[i]*LD.all.11$logFC[i] > 0 & 
        LD.all.12$logFC[i]*LD.all.10$logFC[i] > 0 &
        LD.all.12$logFC[i]*LD.all.9$logFC[i] > 0){
      matched.alt[i]=1
    } else if (LD.all.12$logFC[i]*LD.all.11$logFC[i] > 0 & 
               (LD.all.12$logFC[i]*LD.all.10$logFC[i] > 0 |
                LD.all.12$logFC[i]*LD.all.9$logFC[i] > 0)){
      matched.alt[i]=1
    } else  {
      matched.alt[i]=0
    }
  }else { 
    matched.alt[i]=0
  }
}


trans_temporal.alt= LD.all.12[c(which(matched.alt ==1)),8]
trans_temporal.alt.sig = unique(intersect(trans_temporal.alt,main.LD.all.FDR))
trans_temporal.alt.sig.table = dplyr::filter(LD.all.12, zfin_id_symbol %in% trans_temporal.alt.sig )
trans_temporal.alt.sig.table[,"class"]="trans-temporal"
write.csv(trans_temporal.alt.sig.table, "./outputs/s.table8.trans_temporal_all.LD_DEGs.pn.csv")


#adult_LD_onset, DEG in LD not in (d120 and d13 and d6)

adult_onset.alt = c()

for (i in 1:nrow(LD.all.12)){
  if (abs(LD.all.12$logFC[i])> log2(1.5) & 
      abs(LD.all.10$logFC[i])< log2(1.5) &
      abs(LD.all.9$logFC[i])< log2(1.5)){
    adult_onset.alt[i]=1
  }else { 
    adult_onset.alt[i]=0
  }
}

adult_onset.alt= LD.all.12[c(which(adult_onset.alt == 1)),8]
adult_onset.alt.sig = intersect(adult_onset.alt,LD.all.12$zfin_id_symbol[c(which(LD.all.12$FDR <0.05))])
adult_onset.alt.sig.table = dplyr::filter(LD.all.12, zfin_id_symbol %in% adult_onset.alt.sig )
adult_onset.alt.sig.table[,"class"]="adult-onset"
write.csv(adult_onset.alt.sig.table,"./outputs/s.table8.adult_onset_all.LD_DEGs.pn.csv")

#human_orthologs
combi.LD_deg = c(trans_temporal.alt.sig,adult_onset.alt.sig)
combi.LD_deg.zh = dplyr::filter(zh.all, zfin_id_symbol %in% combi.LD_deg)
#write.csv(combi.dis.LD_deg.zh,"./outputs/combi.dis.LD_deg.zh.csv")

####heatmap trans_temporal.alt.sig

maint.ht.mt = cbind("d6"=LD.all.9$logFC[c(which(LD.all.9$zfin_id_symbol %in% trans_temporal.alt.sig))],
                    "d13"=LD.all.10$logFC[c(which(LD.all.10$zfin_id_symbol %in% trans_temporal.alt.sig))],
                    "d120"=LD.all.11$logFC[c(which(LD.all.11$zfin_id_symbol %in% trans_temporal.alt.sig))],
                    "LD"=LD.all.12$logFC[c(which(LD.all.12$zfin_id_symbol %in% trans_temporal.alt.sig))])
row.names(maint.ht.mt)= LD.all.12$zfin_id_symbol[c(which(LD.all.9$zfin_id_symbol %in% trans_temporal.alt.sig))]

st_sig_tx_1= as.matrix(maint.ht.mt)
distfunc <- function(x) dist(x, method="euclidean")
hclustfunc <- function(x) hclust(x, method="ward.D2")

hc = hclustfunc(distfunc(st_sig_tx_1))
ht.label.order = hc$labels[hc$order]
trans.order=ht.label.order
write(ht.label.order,"./figures/Fig6e_trans_temporal.LD.DEGS2.label.order.pn.txt")
nrow(maint.ht.mt)#479
tiff(paste0("./figures/Fig6e_trans_temporal.LD.DEGS2.heatmap.pn_a",nrow(st_sig_tx_1) ,".tiff"), height = 7,width = 5, res=1200, units = "cm", compression = "lzw")
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

#heatmap for adult_onset.alt.sig

adult_onset.ht.mt = cbind("d6"=LD.all.9$logFC[c(which(LD.all.9$zfin_id_symbol %in% adult_onset.alt.sig))],
                          "d13"=LD.all.10$logFC[c(which(LD.all.10$zfin_id_symbol %in% adult_onset.alt.sig))],
                          "d120"=LD.all.11$logFC[c(which(LD.all.11$zfin_id_symbol %in% adult_onset.alt.sig))],
                          "LD"=LD.all.12$logFC[c(which(LD.all.12$zfin_id_symbol %in% adult_onset.alt.sig))])
row.names(adult_onset.ht.mt)= LD.all.12$zfin_id_symbol[c(which(LD.all.12$zfin_id_symbol %in% adult_onset.alt.sig))]


st_sig_tx_1= as.matrix(adult_onset.ht.mt)
distfunc <- function(x) dist(x, method="euclidean")
hclustfunc <- function(x) hclust(x, method="ward.D2")

hc = hclustfunc(distfunc(st_sig_tx_1))
ht.label.order = hc$labels[hc$order]
adult.order = ht.label.order
write(ht.label.order,"./figures/Fig6e_adult_onset.ELSprimed.LD.DEGS.label.order.pn.txt")
nrow(adult_onset.ht.mt)#1021


tiff(paste0("./figures/Fig6e_adult_onset.ELSprimed.LD.DEGS.2.heatmap.4x4.pn_a",nrow(st_sig_tx_1) ,".tiff"),
     height = 15,width = 5, res=1200, units = "cm", compression = "lzw")

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
              #cellnote=adult_onset.ht.FDR,
              #notecex=1,
              #notecol="black",
              labRow = "",
              keysize = 1 ), silent = T)

dev.off()

#non_transtemportal or adult-onset LD

non_tt.ao = deg245.primed[-c(which(deg245.primed %in% c(adult_onset.alt.sig,trans_temporal.alt.sig)))]

non_tt.ao.LD.ht.mt = cbind("d6"=LD.all.9$logFC[c(which(LD.all.9$zfin_id_symbol %in% non_tt.ao))],
                          "d13"=LD.all.10$logFC[c(which(LD.all.10$zfin_id_symbol %in% non_tt.ao))],
                          "d120"=LD.all.11$logFC[c(which(LD.all.11$zfin_id_symbol %in% non_tt.ao))],
                          "LD"=LD.all.12$logFC[c(which(LD.all.12$zfin_id_symbol %in% non_tt.ao))])
row.names(non_tt.ao.LD.ht.mt)= LD.all.12$zfin_id_symbol[c(which(LD.all.12$zfin_id_symbol %in% non_tt.ao))]


st_sig_tx_1= as.matrix(non_tt.ao.LD.ht.mt)
distfunc <- function(x) dist(x, method="euclidean")
hclustfunc <- function(x) hclust(x, method="ward.D2")

hc = hclustfunc(distfunc(st_sig_tx_1))
ht.label.order = hc$labels[hc$order]
non_tt.ao.order = ht.label.order
write(ht.label.order,"./figures/Fig6e_non_tt.ao.primed.LD.DEGS.label.order.pn.txt")
nrow(non_tt.ao.LD.ht.mt)#747

rhc = hclust(dist(st_sig_tx_1, method="euclidean"), method="ward.D2")
gr.row <- cutree(rhc, 4)
tree.table = data.frame("gene" = names(gr.row), "cluster"= gr.row)
write.csv(tree.table,"./figures/Fig6e.non_tt.ao.primed.LD.DEGS.2.cluster.order.pn.csv")

tiff(paste0("./figures/Fig6e_non_tt.ao.primed.LD.DEGS.2.heatmap.a.pn_",nrow(st_sig_tx_1) ,".tiff"), height = 10,width = 5, res=1200, units = "cm", compression = "lzw")

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
              #RowSideColors=col_col,
              lwid = c(1,2.5),
              #cellnote=adult_onset.ht.FDR,
              #notecex=1,
              #notecol="black",
              labRow = "",
              keysize = 0.5 ), silent = T)

dev.off()

#combined heatmaps

combi_mx = rbind(maint.ht.mt[trans.order,], non_tt.ao.LD.ht.mt[non_tt.ao.order,], adult_onset.ht.mt[adult.order,])
st_sig_tx_1= as.matrix(combi_mx)


tiff(paste0("./figures/Fig6e_non_tt.ao.primed.LD.DEGS.2.heatmap.a.pn_",nrow(st_sig_tx_1) ,".tiff"), height = 6,width = 5, res=1200, units = "cm", compression = "lzw")

try(heatmap.2((st_sig_tx_1), col=mycol, scale='none',
              dendrogram = "none",Colv = F,Rowv = F,
              breaks=breaks,
              density='density', trace = 'none',
              na.color = F,
              offsetRow = 0,
              offsetCol = 0,
              cexCol = 1,
              cexRow = 1,
              #colsep = 4,
              labCol = c("6","13","120","LD"),
              rowsep = c(430,1177),
              sepwidth = c(7,7),
              lwid = c(1,2.5),
              #cellnote=adult_onset.ht.FDR,
              #notecex=1,
              #notecol="black",
              labRow = "",
              keysize = 0.5 ), silent = T)

dev.off()
###combine tables
#gc-prime
LD.deg.list.table[,"GC_primed"]= LD.deg.list.table$zfin_id_symbol
LD.deg.list.table$GC_primed[c(which(LD.deg.list.table$GC_primed %in% deg24.primed))] = "yes"
LD.deg.list.table$GC_primed[c(which(LD.deg.list.table$GC_primed != "yes"))] = "no"

#type
LD.deg.list.table[,"exp.type"]= LD.deg.list.table$zfin_id_symbol

LD.deg.list.table$exp.type[c(which(LD.deg.list.table$exp.type %in% trans_temporal.alt.sig))] = "t.tempo"
LD.deg.list.table$exp.type[c(which(LD.deg.list.table$exp.type %in% adult_onset.alt.sig))] = "ad.onset"

LD.deg.list.table$exp.type[
  intersect(which(LD.deg.list.table$exp.type %in% deg24.primed),
  unique(c(which(LD.deg.list.table$exp.type != "t.tempo"),
                             which(LD.deg.list.table$exp.type !="ad.onset")))
  )
  ] = "non.spec"


write.csv(LD.deg.list.table, "./outputs/s.table3_list of LD-DEGs.pn.csv")


####GO-gprofiler
#read table
trans_go.tb = read.csv("./data/transtemporal_LD_DEG_gProfiler_drerio_05-08-2023_12-47-38__intersections.csv")

#selected term associated genes
#trans-temporal
 top.trans.go.genes = list(
#top BP
"trans.npm"=str_split(trans_go.tb[c(which(trans_go.tb$term_name == "neuron projection morphogenesis")),11],pattern = ",")[[1]],
"trans.csrsp"=str_split(trans_go.tb[c(which(trans_go.tb$term_name == "cell surface receptor signaling pathway")),11],pattern = ",")[[1]],

#topCC
"trans.axon"=str_split(trans_go.tb[c(which(trans_go.tb$term_name == "axon")),11],pattern = ",")[[1]],
"trans.synapse"=str_split(trans_go.tb[c(which(trans_go.tb$term_name == "synapse")),11],pattern = ",")[[1]],

#top MF
"trans.pstk"=str_split(trans_go.tb[c(which(trans_go.tb$term_name == "protein serine/threonine kinase activity")),11],pattern = ",")[[1]]

)

#trans.npm +trans.axon
npm.axon=intersect(top.trans.go.genes[[1]],top.trans.go.genes[[3]])
#trans.npm +trans.synapse
npm.synapse=intersect(top.trans.go.genes[[1]],top.trans.go.genes[[4]])

##trans.csrsp +trans.synapse
csrsp.synapse=intersect(top.trans.go.genes[[2]],top.trans.go.genes[[4]])

##trans.csrsp +trans.pstk
csrsp.pstk=intersect(top.trans.go.genes[[2]],top.trans.go.genes[[5]])
csrsp.npm=intersect(top.trans.go.genes[[2]],top.trans.go.genes[[1]])

deg.cpm12[c(which(deg.cpm12$ensembl_gene_id %in% npm.axon)),c(1,5,8,18,19)]
deg.cpm12[c(which(deg.cpm12$ensembl_gene_id %in% npm.synapse)),c(1,5,8,18,19)]
deg.cpm12[c(which(deg.cpm12$ensembl_gene_id %in% csrsp.synapse)),c(1,5,8,18,19)]
deg.cpm12[c(which(deg.cpm12$ensembl_gene_id %in% csrsp.pstk)),c(1,5,8,18,19)]
deg.cpm12[c(which(deg.cpm12$ensembl_gene_id %in% csrsp.npm)),c(1,5,8,18,19)]



#adult-onset
ad_onset_go.tb = read.csv("./data/adultonset_LD_DEGs_gProfiler_drerio_05-08-2023_12-53-12__intersections.csv")

top.adonset.go.genes = list(
  #top CC
  "adon.memb"=str_split(ad_onset_go.tb[c(which(ad_onset_go.tb$term_name == "membrane")),11],pattern = ",")[[1]]

)

#trans.npm +trans.axon
memb=top.adonset.go.genes[[1]]
deg.cpm12[c(which(deg.cpm12$ensembl_gene_id %in% memb)),c(1,5,8,18,19)]

ad_onset_memb_gene = read.csv("./data/adonset_memb_gProfiler_drerio_20-08-2023_17-14-05__intersections.csv")
memb.genes = list(
  #top BP
  "cell_adhesion"=str_split(ad_onset_memb_gene[c(which(ad_onset_memb_gene$term_name == "cell adhesion")),11],pattern = ",")[[1]],
  "csr_signaling"=str_split(ad_onset_memb_gene[c(which(ad_onset_memb_gene$term_name == "cell surface receptor signaling pathway")),11],pattern = ",")[[1]],
  "synaptic_signaling"=str_split(ad_onset_memb_gene[c(which(ad_onset_memb_gene$term_name == "synaptic signaling")),11],pattern = ",")[[1]]
  
)

#trans.npm +trans.axon
cell_adh=memb.genes[[1]]
csr_signal=memb.genes[[2]]
synaptic_signal=memb.genes[[3]]


deg.cpm12[c(which(deg.cpm12$ensembl_gene_id %in% cell_adh)),c(1,5,8,18,19)]
deg.cpm12[c(which(deg.cpm12$ensembl_gene_id %in% csr_signal)),c(1,5,8,18,19)]
deg.cpm12[c(which(deg.cpm12$ensembl_gene_id %in% synaptic_signal)),c(1,5,8,18,19)]


#### overlapping primed gene bPAC+ vs wt and bPAC+ vs bPAC-


common_primed.deg = intersect(deg12.primed,deg24.primed)
all.primed.deg = unique(c(deg12.primed,deg24.primed))
ft.trans_temporal =  dplyr::filter(trans_temporal.alt.sig.table,zfin_id_symbol %in% common_primed.deg)
ft.adult_onset =  dplyr::filter(adult_onset.alt.sig.table,zfin_id_symbol %in% common_primed.deg)

x = list(
  "bPAC+ vs wt" = deg12.primed,
  "bPAC+ vs bPAC-" = deg24.primed,
  "trans-temporal" = ft.trans_temporal$zfin_id_symbol,
  "adult_onset" = ft.adult_onset$zfin_id_symbol )
 


library(eulerr)

set.seed(3)
fit9= euler(overlapped.primed.list)
fit9

'                                                      original  fitted residuals regionError
bPAC+ vs wt                                                 99  97.923     1.077       0.001
bPAC+ vs bPAC-                                              26   0.000    26.000       0.015
trans-temporal                                               0   3.284    -3.284       0.002
adult_onset                                                  0   3.584    -3.584       0.002
bPAC+ vs wt&bPAC+ vs bPAC-                                 668 670.398    -2.398       0.001
bPAC+ vs wt&trans-temporal                                   0  11.423   -11.423       0.006
bPAC+ vs wt&adult_onset                                      0   9.408    -9.408       0.005
bPAC+ vs bPAC-&trans-temporal                                0   0.000     0.000       0.000
bPAC+ vs bPAC-&adult_onset                                   0   0.000     0.000       0.000
trans-temporal&adult_onset                                   0   0.000     0.000       0.000
bPAC+ vs wt&bPAC+ vs bPAC-&trans-temporal                  312 308.860     3.140       0.003
bPAC+ vs wt&bPAC+ vs bPAC-&adult_onset                     677 675.667     1.333       0.003
bPAC+ vs wt&trans-temporal&adult_onset                       0   0.000     0.000       0.000
bPAC+ vs bPAC-&trans-temporal&adult_onset                    0   0.000     0.000       0.000
bPAC+ vs wt&bPAC+ vs bPAC-&trans-temporal&adult_onset        0  11.206   -11.206       0.006

diagError: 0.015 
stress:    0.001'

#####


# epigenetic regulators ---------------------------------------------------

#downloaded tables from QuickGO with keywords,"epigenetic regulation of gene expression", "DNA modification", "RNA modification","RNA processing", "histone modification"

epi_reg = unique(read.csv("./data/QuickGO_annotations_epigenetic_reg.tsv", sep="\t",header = T)[,3])
DNA_modi = unique(read.csv("./data/QuickGO_annotations_DNA_modification.tsv", sep="\t",header = T)[,3])
hist_modi = unique(read.csv("./data/QuickGO_annotations_histone_modification.tsv", sep="\t",header = T)[,3])
RNA_proc = unique(read.csv("./data/QuickGO_annotations_RNA_processing.tsv", sep="\t",header = T)[,3])
RNA_modi = unique(read.csv("./data/QuickGO_annotations_RNA_modification.tsv", sep="\t",header = T)[,3])


#selection

epi_terms = c("epi_reg","DNA_modi","hist_modi", "RNA_modi","RNA_proc")

for (i in c(9:12)) {
  deg = get(paste0("LD.all.",i))
  deg = deg %>% filter(mean.cpm.sub > 5|mean.cpm.cont >5)%>% filter(FDR < 0.05)%>%
    filter(abs(logFC) >=  log2(1.5))%>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)%>%filter(zfin_id_symbol %in% epi_reg)
  deg[,"epi"] = c(rep(epi_terms[1],nrow(deg)))
  name=paste0("epi_reg_cpm.",i)
  assign(name,deg)
}

for (i in c(9:12)) {
  deg = get(paste0("LD.all.",i))
  deg = deg %>% filter(mean.cpm.sub > 5|mean.cpm.cont >5)%>% filter(FDR < 0.05)%>%
    filter(abs(logFC) >=  log2(1.5))%>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)%>%filter(zfin_id_symbol %in% DNA_modi)
  deg[,"epi"] = c(rep(epi_terms[2],nrow(deg)))
  name=paste0("DNA_modi_cpm.",i)
  assign(name,deg)
}

for (i in c(9:12)) {
  deg = get(paste0("LD.all.",i))
  deg = deg %>% filter(mean.cpm.sub > 5|mean.cpm.cont >5)%>% filter(FDR < 0.05)%>%
    filter(abs(logFC) >=  log2(1.5))%>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)%>%filter(zfin_id_symbol %in% hist_modi)
  deg[,"epi"] = c(rep(epi_terms[3],nrow(deg)))
  name=paste0("hist_modi_cpm.",i)
  assign(name,deg)
}

for (i in c(9:12)) {
  deg = get(paste0("LD.all.",i))
  deg = deg %>% filter(mean.cpm.sub > 5|mean.cpm.cont >5)%>%filter(FDR < 0.05)%>%
    filter(abs(logFC) >=  log2(1.5))%>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)%>%filter(zfin_id_symbol %in% RNA_modi)
  deg[,"epi"] = c(rep(epi_terms[4],nrow(deg)))
  name=paste0("RNA_modi_cpm.",i)
  assign(name,deg)
}


for (i in c(9:12)) {
  deg = get(paste0("LD.all.",i))
  deg = deg %>% filter(mean.cpm.sub > 5|mean.cpm.cont >5)%>% filter(FDR < 0.05)%>%
    filter(abs(logFC) >=  log2(1.5))%>% filter(zfin_id_symbol %in% prot_id_map$SYMBOL)%>%filter(zfin_id_symbol %in% RNA_proc)
  deg[,"epi"] = c(rep(epi_terms[5],nrow(deg)))
  name=paste0("RNA_proc_cpm.",i)
  assign(name,deg)
} 

#trans_temporal GC_affected_sig. epi
m.epi_reg = intersect(unique(c(epi_reg_cpm.9$zfin_id_symbol,epi_reg_cpm.10$zfin_id_symbol)),
                      unique(c(epi_reg_cpm.11$zfin_id_symbol,epi_reg_cpm.12$zfin_id_symbol)))

m.dna.modi = intersect(unique(c(DNA_modi_cpm.9$zfin_id_symbol,DNA_modi_cpm.10$zfin_id_symbol)),
                       unique(c(DNA_modi_cpm.11$zfin_id_symbol,DNA_modi_cpm.12$zfin_id_symbol)))

m.hist_mod = intersect(unique(c(hist_modi_cpm.9$zfin_id_symbol,hist_modi_cpm.10$zfin_id_symbol)),
                       unique(c(hist_modi_cpm.11$zfin_id_symbol,hist_modi_cpm.12$zfin_id_symbol)))

m.RNA.modi = intersect(unique(c(RNA_modi_cpm.9$zfin_id_symbol,RNA_modi_cpm.10$zfin_id_symbol)),
                       unique(c(RNA_modi_cpm.11$zfin_id_symbol,RNA_modi_cpm.12$zfin_id_symbol)))

m.RNA_proc = intersect(unique(c(RNA_proc_cpm.9$zfin_id_symbol,RNA_proc_cpm.10$zfin_id_symbol)),
                       unique(c(RNA_proc_cpm.11$zfin_id_symbol,RNA_proc_cpm.12$zfin_id_symbol)))
#all_sig_early

all_sig_early_epi = rbind(DNA_modi_cpm.9[,-c(subj.cpm.list[[9]],cont.cpm.list[[9]])],
                          DNA_modi_cpm.10[,-c(subj.cpm.list[[10]],cont.cpm.list[[10]])],
                          hist_modi_cpm.9[,-c(subj.cpm.list[[9]],cont.cpm.list[[9]])],
                          hist_modi_cpm.10[,-c(subj.cpm.list[[10]],cont.cpm.list[[10]])],
                          RNA_modi_cpm.9[,-c(subj.cpm.list[[9]],cont.cpm.list[[9]])],
                          RNA_modi_cpm.10[,-c(subj.cpm.list[[10]],cont.cpm.list[[10]])],
                          RNA_proc_cpm.9[,-c(subj.cpm.list[[9]],cont.cpm.list[[9]])],
                          RNA_proc_cpm.10[,-c(subj.cpm.list[[10]],cont.cpm.list[[10]])])


#all_sig_prepostLD

all_sig_ppLD_epi = rbind(DNA_modi_cpm.11[,-c(subj.cpm.list[[11]],cont.cpm.list[[11]])],
                         DNA_modi_cpm.12[,-c(subj.cpm.list[[12]],cont.cpm.list[[12]])],
                         hist_modi_cpm.11[,-c(subj.cpm.list[[11]],cont.cpm.list[[11]])],
                         hist_modi_cpm.12[,-c(subj.cpm.list[[12]],cont.cpm.list[[12]])],
                         RNA_modi_cpm.11[,-c(subj.cpm.list[[11]],cont.cpm.list[[11]])],
                         RNA_modi_cpm.12[,-c(subj.cpm.list[[12]],cont.cpm.list[[12]])],
                         RNA_proc_cpm.11[,-c(subj.cpm.list[[11]],cont.cpm.list[[11]])],
                         RNA_proc_cpm.12[,-c(subj.cpm.list[[12]],cont.cpm.list[[12]])])


LDgene_epi = rbind(all_sig_early_epi,all_sig_ppLD_epi)
write.csv(LDgene_epi,"./outputs/s.table.12_LDgene_epi.pn.csv")

LDgene_epi_dict=c()
LDgene_epi_dict[LDgene_epi$zfin_id_symbol]=LDgene_epi$epi 

  

all_sig_epi = unique(c(all_sig_early_epi$zfin_id_symbol,all_sig_ppLD_epi$zfin_id_symbol))
all_sig_epi_cpm = filter(LD.all.12,zfin_id_symbol %in% all_sig_epi)


#map for epi

all_sig_epi_cpm_table = data.frame()

for(i in c(9:12)){
  deg = get(paste0("deg.cpm",i))
  deg = filter(deg, zfin_id_symbol %in% all_sig_epi)
  deg = deg[,-c(subj.cpm.list[[i]],cont.cpm.list[[i]])]
  all_sig_epi_cpm_table = rbind(deg, all_sig_epi_cpm_table)
}

#matrix for heatmap

sig_epi_ht = as.matrix(data.frame( "d6" = all_sig_epi_cpm_table$logFC[c(which(all_sig_epi_cpm_table$source=="d6pn"))],
                                   "d13"= all_sig_epi_cpm_table$logFC[c(which(all_sig_epi_cpm_table$source=="d13pn"))],
                                   "d120"= all_sig_epi_cpm_table$logFC[c(which(all_sig_epi_cpm_table$source=="d120pn"))],
                                   "LD"= all_sig_epi_cpm_table$logFC[c(which(all_sig_epi_cpm_table$source=="LDpn"))]))
sig_epi_ht.FDR = as.matrix(data.frame( "d6" = all_sig_epi_cpm_table$FDR[c(which(all_sig_epi_cpm_table$source=="d6pn"))],
                                       "d13"= all_sig_epi_cpm_table$FDR[c(which(all_sig_epi_cpm_table$source=="d13pn"))],
                                       "d120"= all_sig_epi_cpm_table$FDR[c(which(all_sig_epi_cpm_table$source=="d120pn"))],
                                       "LD"= all_sig_epi_cpm_table$FDR[c(which(all_sig_epi_cpm_table$source=="LDpn"))]))


sig_epi_ht.FDR[sig_epi_ht.FDR <= 0.05] <-"*"
sig_epi_ht.FDR[sig_epi_ht.FDR > 0.05] <-NA
sig_epi_ht.FDR[abs(sig_epi_ht) < log2(1.5)] <-NA



st_sig_tx_1= sig_epi_ht
distfunc <- function(x) dist(x, method="euclidean")
hclustfunc <- function(x) hclust(x, method="mcquitty")

r.ns = all_sig_epi_cpm_table$zfin_id_symbol[c(which(all_sig_epi_cpm_table$source=="d6pn"))]
row.names(st_sig_tx_1)= r.ns

hc = hclustfunc(distfunc(st_sig_tx_1))
ht.label.order = hc$labels[hc$order]
write(ht.label.order,"./figures/Fig6.epi_LD.label.order.pn.txt")

rhc = hclust(dist(st_sig_tx_1, method="euclidean"), method="mcquitty")
gr.row <- cutree(rhc, 5)
tree.table = data.frame("gene" = names(gr.row), "cluster"= gr.row)
write.csv(tree.table,"./figures/Fig6.epi_LD..cluster.order.pn.csv")
library(viridis)
desturate = viridis(3)
row.col = LDgene_epi_dict[r.ns]
row.col[row.col == "DNA_modi"] <-1
row.col[row.col == "hist_modi"] <-2
row.col[row.col == "RNA_modi"] <-3
row.col[row.col == "RNA_proc"] <-3


tiff(paste0("./figures/Fig6a.epi.LD.DEGS.2.heatmap.pn_",nrow(st_sig_tx_1) ,".tiff"), height = 8,width = 8, res=1200, units = "cm", compression = "lzw")

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


tiff(paste0("./figures/Fig6a.epi.LD.DEGS.2.labeled.heatmap.pn_",nrow(st_sig_tx_1) ,".tiff"), height = 30,width = 10, res=1200, units = "cm", compression = "lzw")

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

epi.bar.table = data.frame("gene" = r.ns, "type" = LDgene_epi_dict[r.ns])

epi.bar= count(epi.bar.table,2)

tiff(paste0("./figures/Fig6.epi.LD.DEGS.2.bar.pn.",nrow(epi.bar.table) ,".tiff"), height = 6,width = 6, res=1200, units = "cm", compression = "lzw")

epi.bar[,"modifications"]= "modulator(53)"
ggplot(epi.bar, aes(x=modifications, y= freq, fill=type)) +
  geom_bar(position="stack",stat="identity") +
  scale_fill_viridis(discrete=TRUE, name="") +
  theme_classic() +
  ylab("") +
  xlab("")
dev.off()



###s.fig7 # comparison with PNAS, https://www.pnas.org/doi/10.1073/pnas.1820842116, Provençal et al., 2020
#Supplementary table 5: List of significantly regulated Pro-diff+WO +acute transcripts that map to a long-lasing DMS (n=702).
s5_pnas = read.csv("./supplementary/pnas.1820842116.sd01_st5.csv")

s5_pnas.fit= na.omit(s5_pnas[,c(5:7)])
s5_pnas.fit=filter(s5_pnas.fit, !grepl(",",s5_pnas.fit$Gene_symbol))
colnames(s5_pnas.fit)=c("human.Symbol","FC","PValue")

homolog_z= filter(zh.all,human.Symbol %in% unique(s5_pnas.fit$human.Symbol))
homolog_z=na.omit(homolog_z)
homolog_z[,"key"]= paste0(homolog_z$human.Symbol,"_",homolog_z$zfin_id_symbol)
d=duplicated(homolog_z$key)
homolog_z=homolog_z[!d,]
homolog_z=na.omit(homolog_z)


library(plyr)


logfc2fc <- function(logFC){
  # sign is -1 if logFC<0; 1 if logFC>=0
  sgn <- (-1)^(1+as.numeric(logFC>=0))
  fc <- sgn*2^abs(logFC)
  return(fc)
}

s5_pnas.fit= na.omit(join(s5_pnas.fit, homolog_z, by="human.Symbol"))


#match
ft.s5.deg246 = deg.cpm12 %>% filter(t.test <0.05) %>% filter(PValue <0.05) %>% filter(mean.cpm.sub > 5|mean.cpm.cont >5) %>%
  filter(zfin_id_symbol %in% homolog_z$zfin_id_symbol)
ft.s5.deg246.fit=ft.s5.deg246[,c(8,1,4)]
ft.s5.deg246.fit[,"FC"]= logfc2fc(ft.s5.deg246.fit$logFC)

ft.s5.deg24 = deg.cpm12 %>% filter(t.test <0.05) %>% filter(PValue <0.05) %>% filter(mean.cpm.sub > 5|mean.cpm.cont >5) %>%
  filter(zfin_id_symbol %in% homolog_z$zfin_id_symbol) %>% filter(zfin_id_symbol %in% deg24.primed)
ft.s5.deg24.fit=ft.s5.deg24[,c(8,1,4)]
ft.s5.deg24.fit[,"FC"]= logfc2fc(ft.s5.deg24.fit$logFC)


ft.s5.deg24.fit=na.omit(join(ft.s5.deg24.fit, homolog_z, by="zfin_id_symbol"))

ft.s5.deg246.fit=na.omit(join(ft.s5.deg246.fit, homolog_z, by="zfin_id_symbol"))

s5_pnas.fit = filter(s5_pnas.fit, key %in% ft.s5.deg246.fit$key)
d=duplicated(s5_pnas.fit$key)
s5_pnas.fit=s5_pnas.fit[!d,]

s5_pnas.fit=arrange(s5_pnas.fit,key)
ft.s5.deg246.fit=arrange(ft.s5.deg246.fit, key)

match= c()
for (i in 1:nrow(ft.s5.deg246.fit)) {
  
  if ((s5_pnas.fit$FC[i]*ft.s5.deg246.fit$FC[i])>0) {
    match=c(match,"overlap")
  }else{
    match = c(match,"no")
  }
}


primed.l= c(ft.s5.deg246.fit$key)
for (i in 1:nrow(ft.s5.deg246.fit)) {
  
  if ((s5_pnas.fit$FC[i]*ft.s5.deg246.fit$FC[i])<0) {
    primed.l[i]=""
  }else if (!(ft.s5.deg246.fit$key[i] %in% ft.s5.deg24.fit$key)){
    primed.l [i]=""
  }
}


df= data.frame("key"=ft.s5.deg246.fit$key,
               "hhc_FC"=s5_pnas.fit$FC,
               "bPAC_FC"= ft.s5.deg246.fit$FC,
               "pval"=sqrt(ft.s5.deg246.fit$PValue*s5_pnas.fit$PValue),
               "match"= match,
               "gc_primed_bPAC"= primed.l)

rownames(df)=ft.s5.deg246.fit$key

library(ggrepel)
set.seed(30)

nw = ggplot(df, aes(hhc_FC,bPAC_FC,label=gc_primed_bPAC,size=-log2(pval)))+
  geom_point(aes(color=match, alpha = 1,size=-log10(pval)))

p1=nw+
  geom_text_repel(max.overlaps = 1000,force = 5,)+
  #geom_text(aes(label = label, size = NULL))+
  theme_classic()+ 
  labs(y="bPAC+LD_FC", x="hpc_ELS-primed_FC")+
  geom_hline(yintercept=c(0,1.13,-1.13),linetype=c("solid",'dotted','dotted'))+
  geom_vline(xintercept=c(0,1.13,-1.13),linetype=c("solid",'dotted','dotted'))+
  #xlim(c(-2.5,2))+ ylim(c(-10,4))+
  scale_color_manual(values = alpha(c("grey","pink"),0.3))+
  theme(axis.title = element_text(size=15, face = "bold"))


tiff(paste0("./figures/s.Fig6.gc-primed.com.pn.tiff"), height = 18,width = 30, res=1200, units = "cm", compression = "lzw")
plot(p1)
dev.off()


write.csv(df,"./outputs/s.table6.GC-prime_gene_comparison.pn.csv")


