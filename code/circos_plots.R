#data
#bPAC+, bPAC-, wt , d6, #bPAC+, bPAC-, wt , d6, d120
#ft.degs

circo.degs.list = list(
"nw.dn"=ft.deg.cpm.dn1$zfin_id_symbol,
"nw.up"=ft.deg.cpm.up1$zfin_id_symbol,

"pw.dn"=ft.deg.cpm.dn5$zfin_id_symbol,
"pw.up"=ft.deg.cpm.up5$zfin_id_symbol,

"np.dn"=ft.deg.cpm.dn9$zfin_id_symbol,
"np.up"=ft.deg.cpm.up9$zfin_id_symbol
)

circo.degs.fc.list = list(
  "nw.dn"=ft.deg.cpm.dn1$logFC,
  "nw.up"=ft.deg.cpm.up1$logFC,
  
  "pw.dn"=ft.deg.cpm.dn5$logFC,
  "pw.up"=ft.deg.cpm.up5$logFC,
  
  "np.dn"=ft.deg.cpm.dn9$logFC,
  "np.up"=ft.deg.cpm.up9$logFC
)



#####data.frame 
samples_df = c("bPAC- vs.wild type.down",
               "bPAC- vs.wild type.up",
           "bPAC+ vs. wild type.down",
           "bPAC+ vs. wild type.up",
           "bAPC+ vs.bPC-.down",
           "bAPC+ vs.bPC-.up")
set.seed(999)
#d6
df_d6 = data.frame("gene_name"=c(circo.degs.list$nw.dn,
                                 rev(circo.degs.list$nw.up),
                                 circo.degs.list$pw.dn,
                                 rev(circo.degs.list$pw.up),
                                 circo.degs.list$np.dn,
                                 rev(circo.degs.list$np.up)), 
                   "sample"=c(rep(samples_df[1],length(circo.degs.list$nw.dn)),
                              rep(samples_df[2],length(circo.degs.list$nw.up)),
                              rep(samples_df[3],length(circo.degs.list$pw.dn)),
                              rep(samples_df[4],length(circo.degs.list$pw.up)),
                              rep(samples_df[5],length(circo.degs.list$np.dn)),
                              rep(samples_df[6],length(circo.degs.list$np.up))), 
                   "logfc"=c(circo.degs.fc.list$nw.dn,
                             rev(circo.degs.fc.list$nw.up),
                             circo.degs.fc.list$pw.dn,
                             rev(circo.degs.fc.list$pw.up),
                             circo.degs.fc.list$np.dn,
                             rev(circo.degs.fc.list$np.up))
                   )

#export
#general track = sector size

deg.num = c(length(circo.degs.list[[1]]),length(circo.degs.list[[2]]),
            length(circo.degs.list[[3]]),length(circo.degs.list[[4]]),
            length(circo.degs.list[[5]]),length(circo.degs.list[[6]]))

track.info = data.frame(
  "sample" = samples_df,
  "start" = 1, 
  "end" = deg.num
)
write.csv(track.info,"./outputs/track.info.csv",row.names = F)

gene_od = c(c(1:deg.num[1]),c(1:deg.num[2]),
            c(1:deg.num[3]),c(1:deg.num[4]),
            c(1:deg.num[5]),c(1:deg.num[6]))

track.logfc = data.frame(
  "sample" = df_d6$sample,
  "start" = gene_od, 
  "end" = gene_od,
  "value1" = df_d6$logfc
)
write.csv(track.logfc,"./outputs/track.logfc.csv",row.names = F)

track.label = data.frame(
  "sample" = df_d6$sample,
  "start" = gene_od, 
  "end" = gene_od,
  "gene_name" = df_d6$gene_name
)
write.csv(track.label,"./outputs/track.label.csv",row.names = F)


  df1 = cbind(df_d6,"start"=gene_od,"end"=gene_od) %>% filter(sample == samples_df[1]) 
  df2 = cbind(df_d6,"start"=gene_od,"end"=gene_od) %>% filter(sample == samples_df[3]) 
  df3 = cbind(df_d6,"start"=gene_od,"end"=gene_od) %>% filter(sample == samples_df[5])
  
  df.int1 = df1 %>% filter(gene_name %in% df2$gene_name)
  df.int2 = df1 %>% filter(gene_name %in% df3$gene_name)
  df.int3 = df2 %>% filter(gene_name %in% df3$gene_name)
  
  df.link1 = left_join(df.int1,df2, by="gene_name")[,c(2,4,5,6,8,9)]
  df.link2 = left_join(df.int2,df3, by="gene_name")[,c(2,4,5,6,8,9)]
  df.link3 = left_join(df.int3,df3, by="gene_name")[,c(2,4,5,6,8,9)]
  df.link1[,"color"]="a"
  df.link2[,"color"]="b"
  df.link3[,"color"]="c"
  
  df4 = cbind(df_d6,"start"=gene_od,"end"=gene_od) %>% filter(sample == samples_df[2]) 
  df5 = cbind(df_d6,"start"=gene_od,"end"=gene_od) %>% filter(sample == samples_df[4]) 
  df6 = cbind(df_d6,"start"=gene_od,"end"=gene_od) %>% filter(sample == samples_df[6])
  
  df.int4 = df4 %>% filter(gene_name %in% df5$gene_name)
  df.int5 = df4 %>% filter(gene_name %in% df6$gene_name)
  df.int6 = df5 %>% filter(gene_name %in% df6$gene_name)
  
  df.link4 = left_join(df.int4,df5, by="gene_name")[,c(2,4,5,6,8,9)]
  df.link5 = left_join(df.int5,df6, by="gene_name")[,c(2,4,5,6,8,9)]
  df.link6 = left_join(df.int6,df6, by="gene_name")[,c(2,4,5,6,8,9)]
  df.link4[,"color"]="d"
  df.link5[,"color"]="e"
  df.link6[,"color"]="f"
  
  

track.link1 = rbind(df.link1,df.link2,df.link3)
track.link2 = rbind(df.link4,df.link5,df.link6)

write.csv(track.link1,"./outputs/track.link1.csv",row.names = F)
write.csv(track.link2,"./outputs/track.link2.csv",row.names = F)

viridis(6)
#a:#440154FF;b:#414487FF;c:#2A788EFF;d:#22A884FF;e:#7AD151FF;f:#FDE725FF

##440154FF,grey,#414487FF,grey,#2A788EFF,grey

#grey,#414487FF,grey,#22A884FF,grey,#FDE725FF
####
#library(shiny)
#runApp("d:/projects/shinyCircos-master/shinyCircos-master/", launch.browser = TRUE)




#data
#bPAC+, bPAC-, wt , d6, #bPAC+, bPAC-, wt , d6, d120
#ft.degs

circo.degs.list = list(
  "nw.dn"=ft.deg.cpm.dn2$zfin_id_symbol,
  "nw.up"=ft.deg.cpm.up2$zfin_id_symbol,
  
  "pw.dn"=ft.deg.cpm.dn6$zfin_id_symbol,
  "pw.up"=ft.deg.cpm.up6$zfin_id_symbol,
  
  "np.dn"=ft.deg.cpm.dn10$zfin_id_symbol,
  "np.up"=ft.deg.cpm.up10$zfin_id_symbol
)

circo.degs.fc.list = list(
  "nw.dn"=ft.deg.cpm.dn2$logFC,
  "nw.up"=ft.deg.cpm.up2$logFC,
  
  "pw.dn"=ft.deg.cpm.dn6$logFC,
  "pw.up"=ft.deg.cpm.up6$logFC,
  
  "np.dn"=ft.deg.cpm.dn10$logFC,
  "np.up"=ft.deg.cpm.up10$logFC
)



#####data.frame 
samples_df = c("bPAC- vs.wild type.down",
               "bPAC- vs.wild type.up",
               "bPAC+ vs. wild type.down",
               "bPAC+ vs. wild type.up",
               "bAPC+ vs.bPC-.down",
               "bAPC+ vs.bPC-.up")
set.seed(999)
#d13
df_d13 = data.frame("gene_name"=c(circo.degs.list$nw.dn,
                                  rev(circo.degs.list$nw.up),
                                  circo.degs.list$pw.dn,
                                  rev(circo.degs.list$pw.up),
                                  circo.degs.list$np.dn,
                                  rev(circo.degs.list$np.up)), 
                    "sample"=c(rep(samples_df[1],length(circo.degs.list$nw.dn)),
                               rep(samples_df[2],length(circo.degs.list$nw.up)),
                               rep(samples_df[3],length(circo.degs.list$pw.dn)),
                               rep(samples_df[4],length(circo.degs.list$pw.up)),
                               rep(samples_df[5],length(circo.degs.list$np.dn)),
                               rep(samples_df[6],length(circo.degs.list$np.up))), 
                    "logfc"=c(circo.degs.fc.list$nw.dn,
                              rev(circo.degs.fc.list$nw.up),
                              circo.degs.fc.list$pw.dn,
                              rev(circo.degs.fc.list$pw.up),
                              circo.degs.fc.list$np.dn,
                              rev(circo.degs.fc.list$np.up))
)

#export
#general track = sector size

deg.num = c(length(circo.degs.list[[1]]),length(circo.degs.list[[2]]),
            length(circo.degs.list[[3]]),length(circo.degs.list[[4]]),
            length(circo.degs.list[[5]]),length(circo.degs.list[[6]]))

track.info = data.frame(
  "sample" = samples_df,
  "start" = 1, 
  "end" = deg.num
)
write.csv(track.info,"./outputs/track.info13.csv",row.names = F)

gene_od = c(c(1:deg.num[1]),c(1:deg.num[2]),
            c(1:deg.num[3]),c(1:deg.num[4]),
            c(1:deg.num[5]),c(1:deg.num[6]))

track.logfc = data.frame(
  "sample" = df_d13$sample,
  "start" = gene_od, 
  "end" = gene_od,
  "value1" = df_d13$logfc
)
write.csv(track.logfc,"./outputs/track.logfc13.csv",row.names = F)

track.label = data.frame(
  "sample" = df_d13$sample,
  "start" = gene_od, 
  "end" = gene_od,
  "gene_name" = df_d13$gene_name
)
write.csv(track.label,"./outputs/track.label13.csv",row.names = F)


df1 = cbind(df_d13,"start"=gene_od,"end"=gene_od) %>% filter(sample == samples_df[1]) 
df2 = cbind(df_d13,"start"=gene_od,"end"=gene_od) %>% filter(sample == samples_df[3]) 
df3 = cbind(df_d13,"start"=gene_od,"end"=gene_od) %>% filter(sample == samples_df[5])

df.int1 = df1 %>% filter(gene_name %in% df2$gene_name)
df.int2 = df1 %>% filter(gene_name %in% df3$gene_name)
df.int3 = df2 %>% filter(gene_name %in% df3$gene_name)

df.link1 = left_join(df.int1,df2, by="gene_name")[,c(2,4,5,6,8,9)]
df.link2 = left_join(df.int2,df3, by="gene_name")[,c(2,4,5,6,8,9)]
df.link3 = left_join(df.int3,df3, by="gene_name")[,c(2,4,5,6,8,9)]
df.link1[,"color"]="a"
df.link2[,"color"]="b"
df.link3[,"color"]="c"

df4 = cbind(df_d13,"start"=gene_od,"end"=gene_od) %>% filter(sample == samples_df[2]) 
df5 = cbind(df_d13,"start"=gene_od,"end"=gene_od) %>% filter(sample == samples_df[4]) 
df6 = cbind(df_d13,"start"=gene_od,"end"=gene_od) %>% filter(sample == samples_df[6])

df.int4 = df4 %>% filter(gene_name %in% df5$gene_name)
df.int5 = df4 %>% filter(gene_name %in% df6$gene_name)
df.int6 = df5 %>% filter(gene_name %in% df6$gene_name)

df.link4 = left_join(df.int4,df5, by="gene_name")[,c(2,4,5,6,8,9)]
df.link5 = left_join(df.int5,df6, by="gene_name")[,c(2,4,5,6,8,9)]
df.link6 = left_join(df.int6,df6, by="gene_name")[,c(2,4,5,6,8,9)]
df.link4[,"color"]="d"
df.link5[,"color"]="e"
df.link6[,"color"]="f"



track.link1 = rbind(df.link1,df.link2,df.link3)
track.link2 = rbind(df.link4,df.link5,df.link6)

write.csv(track.link1,"./outputs/track.link1_13.csv",row.names = F)
write.csv(track.link2,"./outputs/track.link2_13.csv",row.names = F)

viridis(6)

#a:#440154FF;b:#414487FF;c:#2A788EFF;d:#22A884FF;e:#7AD151FF;f:#FDE725FF

##440154FF,grey,#414487FF,grey,#2A788EFF,grey

#grey,#414487FF,grey,#22A884FF,#grey,#FDE725FF

####
#library(shiny)
#runApp("d:/projects/shinyCircos-master/shinyCircos-master/", launch.browser = TRUE)




#data
#bPAC+, bPAC-, wt , d6, #bPAC+, bPAC-, wt , d6, d120
#ft.degs

circo.degs.list = list(
  "nw.dn"=ft.deg.cpm.dn3$zfin_id_symbol,
  "nw.up"=ft.deg.cpm.up3$zfin_id_symbol,
  
  "pw.dn"=ft.deg.cpm.dn7$zfin_id_symbol,
  "pw.up"=ft.deg.cpm.up7$zfin_id_symbol,
  
  "np.dn"=ft.deg.cpm.dn11$zfin_id_symbol,
  "np.up"=ft.deg.cpm.up11$zfin_id_symbol
)

circo.degs.fc.list = list(
  "nw.dn"=ft.deg.cpm.dn3$logFC,
  "nw.up"=ft.deg.cpm.up3$logFC,
  
  "pw.dn"=ft.deg.cpm.dn7$logFC,
  "pw.up"=ft.deg.cpm.up7$logFC,
  
  "np.dn"=ft.deg.cpm.dn11$logFC,
  "np.up"=ft.deg.cpm.up11$logFC
)



#####data.frame 
samples_df = c("bPAC- vs.wild type.down",
               "bPAC- vs.wild type.up",
               "bPAC+ vs. wild type.down",
               "bPAC+ vs. wild type.up",
               "bAPC+ vs.bPC-.down",
               "bAPC+ vs.bPC-.up")
set.seed(999)
#d120
df_d120 = data.frame("gene_name"=c(circo.degs.list$nw.dn,
                                   rev(circo.degs.list$nw.up),
                                   circo.degs.list$pw.dn,
                                   rev(circo.degs.list$pw.up),
                                   circo.degs.list$np.dn,
                                   rev(circo.degs.list$np.up)), 
                     "sample"=c(rep(samples_df[1],length(circo.degs.list$nw.dn)),
                                rep(samples_df[2],length(circo.degs.list$nw.up)),
                                rep(samples_df[3],length(circo.degs.list$pw.dn)),
                                rep(samples_df[4],length(circo.degs.list$pw.up)),
                                rep(samples_df[5],length(circo.degs.list$np.dn)),
                                rep(samples_df[6],length(circo.degs.list$np.up))), 
                     "logfc"=c(circo.degs.fc.list$nw.dn,
                               rev(circo.degs.fc.list$nw.up),
                               circo.degs.fc.list$pw.dn,
                               rev(circo.degs.fc.list$pw.up),
                               circo.degs.fc.list$np.dn,
                               rev(circo.degs.fc.list$np.up))
)

#export
#general track = sector size

deg.num = c(length(circo.degs.list[[1]]),length(circo.degs.list[[2]]),
            length(circo.degs.list[[3]]),length(circo.degs.list[[4]]),
            length(circo.degs.list[[5]]),length(circo.degs.list[[6]]))

track.info = data.frame(
  "sample" = samples_df,
  "start" = 1, 
  "end" = deg.num
)
write.csv(track.info,"./outputs/track.info120.csv",row.names = F)

gene_od = c(c(1:deg.num[1]),c(1:deg.num[2]),
            c(1:deg.num[3]),c(1:deg.num[4]),
            c(1:deg.num[5]),c(1:deg.num[6]))

track.logfc = data.frame(
  "sample" = df_d120$sample,
  "start" = gene_od, 
  "end" = gene_od,
  "value1" = df_d120$logfc
)
write.csv(track.logfc,"./outputs/track.logfc120.csv",row.names = F)

track.label = data.frame(
  "sample" = df_d120$sample,
  "start" = gene_od, 
  "end" = gene_od,
  "gene_name" = df_d120$gene_name
)
write.csv(track.label,"./outputs/track.label120.csv",row.names = F)


df1 = cbind(df_d120,"start"=gene_od,"end"=gene_od) %>% filter(sample == samples_df[1]) 
df2 = cbind(df_d120,"start"=gene_od,"end"=gene_od) %>% filter(sample == samples_df[3]) 
df3 = cbind(df_d120,"start"=gene_od,"end"=gene_od) %>% filter(sample == samples_df[5])

df.int1 = df1 %>% filter(gene_name %in% df2$gene_name)
df.int2 = df1 %>% filter(gene_name %in% df3$gene_name)
df.int3 = df2 %>% filter(gene_name %in% df3$gene_name)

df.link1 = left_join(df.int1,df2, by="gene_name")[,c(2,4,5,6,8,9)]
df.link2 = left_join(df.int2,df3, by="gene_name")[,c(2,4,5,6,8,9)]
df.link3 = left_join(df.int3,df3, by="gene_name")[,c(2,4,5,6,8,9)]
df.link1[,"color"]="a"
df.link2[,"color"]="b"
df.link3[,"color"]="c"

df4 = cbind(df_d120,"start"=gene_od,"end"=gene_od) %>% filter(sample == samples_df[2]) 
df5 = cbind(df_d120,"start"=gene_od,"end"=gene_od) %>% filter(sample == samples_df[4]) 
df6 = cbind(df_d120,"start"=gene_od,"end"=gene_od) %>% filter(sample == samples_df[6])

df.int4 = df4 %>% filter(gene_name %in% df5$gene_name)
df.int5 = df4 %>% filter(gene_name %in% df6$gene_name)
df.int6 = df5 %>% filter(gene_name %in% df6$gene_name)

df.link4 = left_join(df.int4,df5, by="gene_name")[,c(2,4,5,6,8,9)]
df.link5 = left_join(df.int5,df6, by="gene_name")[,c(2,4,5,6,8,9)]
df.link6 = left_join(df.int6,df6, by="gene_name")[,c(2,4,5,6,8,9)]
df.link4[,"color"]="d"
df.link5[,"color"]="e"
df.link6[,"color"]="f"



track.link1 = rbind(df.link1,df.link2,df.link3)
track.link2 = rbind(df.link4,df.link5,df.link6)

write.csv(track.link1,"./outputs/track.link1_120.csv",row.names = F)
write.csv(track.link2,"./outputs/track.link2_120.csv",row.names = F)

viridis(6)
#a:#440154FF;b:#414487FF;c:#2A788EFF;d:#22A884FF;e:#7AD151FF;f:#FDE725FF

##440154FF,grey,#414487FF,grey,#2A788EFF,grey

#grey,#414487FF,grey,#22A884FF,#grey,#FDE725FF
####
#library(shiny)
#runApp("d:/projects/shinyCircos-master/shinyCircos-master/", launch.browser = TRUE)