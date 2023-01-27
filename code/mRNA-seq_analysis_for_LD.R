#DESeq2_DEG_analysis

#countr.matrix.LD for all

#countr.matrix.LD for TU and STAR

countr.matrix_LD = makeContrasts(
  
  d120_s.POSvsWT = star_pos_d120 - tu_wt_d120,
  d120_LD_s.POSvsWT = star_pos_LD - tu_wt_LD,
  d120_s.POSvsd120_LD_s.POS = star_pos_LD - star_pos_d120,
  d120_TUvsd120_LD.TU = tu_wt_LD - tu_wt_d120,
  
  levels = colnames(design))
countr.matrix_LD



##Removing heteroscedascity
#
gfit.LD <- glmQLFit(ft_x.2.norm, design, robust = T)
c.gfit.LD = glmQLFTest(gfit, contrast=countr.matrix_LD)

###3
call.name = colnames(countr.matrix_LD)
for(i in c(1:4)){
  qn = paste0("qlf.LD",i)  
  assign(qn, glmQLFTest(gfit, contrast=countr.matrix_LD[,call.name[i]]))
  
  print(i)
}

#test
for(i in c(1:4)){
  qlf = get(paste0("qlf.LD",i))
  FDR = paste0("FDR.LD",i)
  assign(FDR, p.adjust(qlf$table$PValue, method="BH"))
  
  print(i)
}



gfit.ens = ID.gene[c(which(ID.gene$id %in% rownames(gfit.LD))),]
fit.norm.counts = gfit.LD$counts
all_count = cbind(gfit.ens,fit.norm.counts)
cmp_all = cpm(ft_x.2.norm)


for(i in 1:4){
  qn = paste0("qlf.LD",i)
  qt = topTags(get(qn),n = Inf)
  assign(paste0("topqlf.LD",i), data.frame(qt$table))
}

for(i in 1:4){
  qn = get(paste0("topqlf.LD",i))
  qn = qn[gfit.ens$id,]
  assign(paste0("topqlf.LD",i), qn)
}

#d120_s.POSvsd120_LD_s.POS=star_pos_d120_LD-star_pos_d120,
deg.cpm9= cbind(topqlf.LD3,gfit.ens,cmp_all[,c(21:30)])
write.csv(deg.cpm9,"./data/topqlf.cpm.LD.3.csv")


#d120_TUvsd120_LD.TU=tu_wt_d120_LD-tu_wt_d120,
deg.cpm10= cbind(topqlf.LD4,gfit.ens,cmp_all[,c(41:50)])
write.csv(deg.cpm10,"./data/topqlf.cpm.LD.4.csv")






#fc+counts
#1:d120_s.POSvsWT = star_pos_d120 - tu_wt_d120,
#2:d120_LD_s.POSvsWT = star_pos_d120_LD - tu_wt_d120_LD,
#3:d120_s.POSvsd120_LD_s.POS = star_pos_d120_LD - star_pos_d120,
#4:d120_TUvsd120_LD.TU = tu_wt_d120_LD - tu_wt_d120,

