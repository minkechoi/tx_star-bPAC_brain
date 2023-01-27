
# GO_analyses -------------------------------------------------------------

#list of background genes
for (i in c(1:3,5:7)){
  deg = get(paste0("deg.cpm",i))
  deg = deg %>% filter(mean.cpm.pos > 5|mean.cpm.wt >5)
  deg = prot_id_map[deg$zfin_id_symbol, "ENTREZID"]
  deg = na.omit(deg)
  name = paste0("bk.gene",i)
  assign(name, deg)
  print(i)
}

bk.all.basal = unique(c(bk.gene1,bk.gene2,bk.gene3,bk.gene5,bk.gene6,bk.gene7))
###bPAC+ time course 

#up degs
tc_point_up <- list()
## Reformat the timecourse cluster classifications as a list

#altering gene ID to ENTREZID

for (i in c(5:7)){
  cluster_psites <- which(deg_all.up$source == time.name[i])
  cluster_prots <- unique(deg_all.up$zfin_id_symbol[cluster_psites])
  cluster_prots_entrez <- prot_id_map[cluster_prots, "ENTREZID"]
  cluster_prots_entrez=na.omit(cluster_prots_entrez)
  tc_point_up[[time.name[i]]] <- cluster_prots_entrez
}
## Do the enrichment analysis using clusterProfiler
##GOBP
tc_enrich.up <- clusterProfiler::compareCluster(tc_point_up, fun = "enrichGO",  OrgDb='org.Dr.eg.db', ont = "BP",
                                                bk.all.basal,
                                                pvalueCutoff = 0.05,
                                                pAdjustMethod = "BH",qvalueCutoff = 0.2,
                                                minGSSize = 10,
                                                maxGSSize = 600, 
                                                readable = TRUE)

## tables
tc_enrich.up.bp.table = data.frame(tc_enrich.up[1:length(tc_enrich.up[]$ID)])
#write.csv(tc_enrich.up.bp.table,"./outputs/tc_enrich.up.bp.table.csv")

####GOMF
tc_enrich.up.mf <- clusterProfiler::compareCluster(tc_point_up, fun = "enrichGO",  OrgDb='org.Dr.eg.db', ont = "MF",
                                                   bk.all.basal,
                                                   pvalueCutoff = 0.05,
                                                   pAdjustMethod = "BH",qvalueCutoff = 0.2,
                                                   minGSSize = 10,
                                                   maxGSSize = 600, 
                                                   readable = TRUE)

## tables
tc_enrich.up.mf.table = data.frame(tc_enrich.up.mf[1:length(tc_enrich.up.mf[]$ID)])
#write.csv(tc_enrich.up.mf.table,"./outputs/tc_enrich.up.mf.table.csv")

###GOCC
tc_enrich.up.cc <- clusterProfiler::compareCluster(tc_point_up, fun = "enrichGO",  OrgDb='org.Dr.eg.db', ont = "CC",
                                                   bk.all.basal,pvalueCutoff = 0.05,
                                                   pAdjustMethod = "BH",qvalueCutoff = 0.2,
                                                   minGSSize = 10,
                                                   maxGSSize = 600, 
                                                   readable = TRUE)

## tables
tc_enrich.up.cc.table = data.frame(tc_enrich.up.cc[1:length(tc_enrich.up.cc[]$ID)])
#write.csv(tc_enrich.up.cc.table,"./outputs/tc_enrich.up.cc.table.csv")


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
                                                maxGSSize = 1000, readable = TRUE)

##tables

tc_enrich.dn.bp.table = data.frame(tc_enrich.dn[1:length(tc_enrich.dn[]$ID)])
#write.csv(tc_enrich.dn.bp.table,"./outputs/tc_enrich.dn.bp.table.csv")

#write.csv(tc_enrich.dn.table,"./tc_enrich.dn.table.csv")
# can adjust the number of representing terms, showCategory=5

####GOMF
tc_enrich.dn.mf <- clusterProfiler::compareCluster(tc_point_dn, fun = "enrichGO",  OrgDb='org.Dr.eg.db', ont = "MF",
                                                   bk.all.basal,pvalueCutoff = 0.05,
                                                   pAdjustMethod = "BH",qvalueCutoff = 0.2,
                                                   minGSSize = 10,
                                                   maxGSSize = 1000, 
                                                   readable = TRUE)


##tables

tc_enrich.dn.mf.table = data.frame(tc_enrich.dn.mf[1:length(tc_enrich.dn.mf[]$ID)])
#write.csv(tc_enrich.dn.mf.table,"./outputs/tc_enrich.dn.mf.table.csv")

###GOCC
tc_enrich.dn.cc <- clusterProfiler::compareCluster(tc_point_dn, fun = "enrichGO",  OrgDb='org.Dr.eg.db', ont = "CC",
                                                   bk.all.basal,pvalueCutoff = 0.05,
                                                   pAdjustMethod = "BH",qvalueCutoff = 0.2,
                                                   minGSSize = 10,
                                                   maxGSSize = 1000, 
                                                   readable = TRUE)


##tables
tc_enrich.dn.cc.table = data.frame(tc_enrich.dn.cc[1:length(tc_enrich.dn.cc[]$ID)])
#write.csv(tc_enrich.dn.cc.table,"./outputs/tc_enrich.dn.cc.table.csv")

#generating tables

gos=c("bp","mf","cc")
for (i in c(1:3)) {
  #for up-reg.
  tb.up=get(paste0("tc_enrich.up.",gos[i],".table"))
  tb.up[,"source"]=paste0(gos[i],"_up")
 
  a= str_split(tb.up[,4],pattern = "/")
  if (length(a)>1) {
    b= data.frame(Reduce(rbind,a))
  } else {b= data.frame(a[[1]][1],a[[1]][2])  
  }  
  tb.up[,"n.deg"] = as.numeric(b[,2])
 
  a= str_split(tb.up[,5],pattern = "/")
  if (length(a)>1) {
    b= data.frame(Reduce(rbind,a))
  } else {b= data.frame(a[[1]][1],a[[1]][2])  
  }  
  tb.up[,"gene.in.catg"] = as.numeric(b[,1])
  
  tb.up[,"ratio.catg"] = tb.up[,10]/tb.up[,"gene.in.catg"]
  
  #filtering
  tb.up= tb.up %>% filter(Count >= 10) %>% filter(ratio.catg >= 0.09)
  
  
  #for down-reg  
  tb.dn=get(paste0("tc_enrich.dn.",gos[i],".table"))
  tb.dn[,"source"]=paste0(gos[i],"_dn")

  a= str_split(tb.dn[,4],pattern = "/")
  if (length(a)>1) {
    b= data.frame(Reduce(rbind,a))
  } else {b= data.frame(a[[1]][1],a[[1]][2])  
  }  
  tb.dn[,"n.deg"] = as.numeric(b[,2])
  
  a= str_split(tb.dn[,5],pattern = "/")
  if (length(a)>1) {
    b= data.frame(Reduce(rbind,a))
  } else {b= data.frame(a[[1]][1],a[[1]][2])  
  }  
  tb.dn[,"gene.in.catg"] = as.numeric(b[,1])
  
  tb.dn[,"ratio.catg"] = c(tb.dn[,10]/tb.dn[,"gene.in.catg"])
  
  #filtering
  tb.dn= tb.dn %>% filter(Count >= 10) %>% filter(ratio.catg >= 0.09)
  
  ###
    
  b.tb=rbind(tb.up,tb.dn)
  name= paste0("tc_enrich_sp_",gos[i])
  assign(name ,b.tb )
}

#write.csv(tc_enrich_sp_bp,"./tc_enrich_sp_bp.csv")
#write.csv(tc_enrich_sp_mf,"./tc_enrich_sp_mf.csv")
#write.csv(tc_enrich_sp_cc,"./tc_enrich_sp_cc.csv")
tc_enrich_sp= rbind(tc_enrich_sp_bp,tc_enrich_sp_mf,tc_enrich_sp_cc)
write.csv(tc_enrich_sp,"./outputs/s.table2a_tc_enrich_sp.csv")

#selection of GOterms



###bPAC- and + individual timepoint analysis

#up degs
tc_point_up <- list()
## Reformat the timecourse cluster classifications as a list

#altering gene ID to ENTREZID

for (i in c(1:3)){
  cluster_psites <- which(deg_all.up$source == time.name[i])
  cluster_prots <- unique(deg_all.up$zfin_id_symbol[cluster_psites])
  cluster_prots_entrez <- prot_id_map[cluster_prots, "ENTREZID"]
  cluster_prots_entrez=na.omit(cluster_prots_entrez)
  tc_point_up[[time.name[i]]] <- cluster_prots_entrez
}
## Do the enrichment analysis using clusterProfiler
##GOBP
tc_enrich.up <- clusterProfiler::compareCluster(tc_point_up, fun = "enrichGO",  OrgDb='org.Dr.eg.db', ont = "BP",
                                                bk.all.basal,
                                                pvalueCutoff = 0.05,
                                                pAdjustMethod = "BH",qvalueCutoff = 0.2,
                                                minGSSize = 10,
                                                maxGSSize = 600, 
                                                readable = TRUE)

## tables
tc_enrich.up.bp.table = data.frame(tc_enrich.up[1:length(tc_enrich.up[]$ID)])
#write.csv(tc_enrich.up.bp.table,"./outputs/tc_enrich.up.bp.table.csv")

####GOMF
tc_enrich.up.mf <- clusterProfiler::compareCluster(tc_point_up, fun = "enrichGO",  OrgDb='org.Dr.eg.db', ont = "MF",
                                                   bk.all.basal,
                                                   pvalueCutoff = 0.05,
                                                   pAdjustMethod = "BH",qvalueCutoff = 0.2,
                                                   minGSSize = 10,
                                                   maxGSSize = 600, 
                                                   readable = TRUE)

## tables
tc_enrich.up.mf.table = data.frame(tc_enrich.up.mf[1:length(tc_enrich.up.mf[]$ID)])
#write.csv(tc_enrich.up.mf.table,"./outputs/tc_enrich.up.mf.table.csv")

###GOCC
tc_enrich.up.cc <- clusterProfiler::compareCluster(tc_point_up, fun = "enrichGO",  OrgDb='org.Dr.eg.db', ont = "CC",
                                                   bk.all.basal,pvalueCutoff = 0.05,
                                                   pAdjustMethod = "BH",qvalueCutoff = 0.2,
                                                   minGSSize = 10,
                                                   maxGSSize = 600, 
                                                   readable = TRUE)

## tables
tc_enrich.up.cc.table = data.frame(tc_enrich.up.cc[1:length(tc_enrich.up.cc[]$ID)])
#write.csv(tc_enrich.up.cc.table,"./outputs/tc_enrich.up.cc.table.csv")


#down degs
tc_point_dn <- list()
## Reformat the timecourse cluster classifications as a list

for (i in c(1:3)){
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
                                                maxGSSize = 1000, readable = TRUE)

##tables

tc_enrich.dn.bp.table = data.frame(tc_enrich.dn[1:length(tc_enrich.dn[]$ID)])
#write.csv(tc_enrich.dn.bp.table,"./outputs/tc_enrich.dn.bp.table.csv")

#write.csv(tc_enrich.dn.table,"./tc_enrich.dn.table.csv")
# can adjust the number of representing terms, showCategory=5

####GOMF
tc_enrich.dn.mf <- clusterProfiler::compareCluster(tc_point_dn, fun = "enrichGO",  OrgDb='org.Dr.eg.db', ont = "MF",
                                                   bk.all.basal,pvalueCutoff = 0.05,
                                                   pAdjustMethod = "BH",qvalueCutoff = 0.2,
                                                   minGSSize = 10,
                                                   maxGSSize = 1000, 
                                                   readable = TRUE)


##tables

tc_enrich.dn.mf.table = data.frame(tc_enrich.dn.mf[1:length(tc_enrich.dn.mf[]$ID)])
#write.csv(tc_enrich.dn.mf.table,"./outputs/tc_enrich.dn.mf.table.csv")

###GOCC
tc_enrich.dn.cc <- clusterProfiler::compareCluster(tc_point_dn, fun = "enrichGO",  OrgDb='org.Dr.eg.db', ont = "CC",
                                                   bk.all.basal,pvalueCutoff = 0.05,
                                                   pAdjustMethod = "BH",qvalueCutoff = 0.2,
                                                   minGSSize = 10,
                                                   maxGSSize = 1000, 
                                                   readable = TRUE)


##tables
tc_enrich.dn.cc.table = data.frame(tc_enrich.dn.cc[1:length(tc_enrich.dn.cc[]$ID)])
#write.csv(tc_enrich.dn.cc.table,"./outputs/tc_enrich.dn.cc.table.csv")

#generating tables

gos=c("bp","mf","cc")
for (i in c(1:3)) {
  #for up-reg.
  tb.up=get(paste0("tc_enrich.up.",gos[i],".table"))
  tb.up[,"source"]=paste0(gos[i],"_up")
  
  a= str_split(tb.up[,4],pattern = "/")
  if (length(a)>1) {
    b= data.frame(Reduce(rbind,a))
  } else {b= data.frame(a[[1]][1],a[[1]][2])  
  }  
  tb.up[,"n.deg"] = as.numeric(b[,2])
  
  a= str_split(tb.up[,5],pattern = "/")
  if (length(a)>1) {
    b= data.frame(Reduce(rbind,a))
  } else {b= data.frame(a[[1]][1],a[[1]][2])  
  }  
  tb.up[,"gene.in.catg"] = as.numeric(b[,1])
  
  tb.up[,"ratio.catg"] = tb.up[,10]/tb.up[,"gene.in.catg"]
  
  #filtering
  #tb.up= tb.up %>% filter(Count >= 10) %>% filter(ratio.catg >= 0.09)
  
  
  #for down-reg  
  tb.dn=get(paste0("tc_enrich.dn.",gos[i],".table"))
  tb.dn[,"source"]=paste0(gos[i],"_dn")
  
  a= str_split(tb.dn[,4],pattern = "/")
  if (length(a)>1) {
    b= data.frame(Reduce(rbind,a))
  } else {b= data.frame(a[[1]][1],a[[1]][2])  
  }  
  tb.dn[,"n.deg"] = as.numeric(b[,2])
  
  a= str_split(tb.dn[,5],pattern = "/")
  if (length(a)>1) {
    b= data.frame(Reduce(rbind,a))
  } else {b= data.frame(a[[1]][1],a[[1]][2])  
  }  
  tb.dn[,"gene.in.catg"] = as.numeric(b[,1])
  
  tb.dn[,"ratio.catg"] = c(tb.dn[,10]/tb.dn[,"gene.in.catg"])
  
  #filtering
  #tb.dn= tb.dn %>% filter(Count >= 3) %>% filter(ratio.catg >= 0.01)
  
  ###
  
  b.tb=rbind(tb.up,tb.dn)
  name= paste0("tc_enrich_sn_",gos[i])
  assign(name ,b.tb )
}

#write.csv(tc_enrich_sn_bp,"./tc_enrich_sn_bp.csv")
#write.csv(tc_enrich_sn_mf,"./tc_enrich_sn_mf.csv")
#write.csv(tc_enrich_sn_cc,"./tc_enrich_sn_cc.csv")
tc_enrich_sn= rbind(tc_enrich_sn_bp,tc_enrich_sn_mf,tc_enrich_sn_cc)
write.csv(tc_enrich_sn,"./outputs/s.table2b_tc_enrich_sn.csv")




