
##read homologs of human and mice ----------------------------------------


zh.all = read.csv("./data/zh.all.csv",header = T)
zm.all = read.csv("./data/zm.all.csv",header = T)


##known oxt regulatory elements ------------------------------------------


oxt_reg = c("oxt","oxtr","otpa","otpb","myt1la","fez2","sim1a","arnt2","olig2","pou3f2a","med12","nr2f2","lhx5","lef1","fezf2","esrrb","thrb","thraa","thrab","rargb","raraa","rarab","rarga","nr2c2","nr2f1b","nr2f1a")



# ##### genes list from QiuckGO: keyword: social behavior -----------------

sb_mus = read.table("./data/social_behavior_mus.tsv",sep = "\t",header = T)
sb_homo = read.table("./data/social_behavior_homo.tsv",sep = "\t",header = T)


# #identifying danio homologs ---------------------------------------------


sb_zm = filter(zm.all, Mouse.Symbol %in% sb_mus$SYMBOL)  
sb_zh = filter(zh.all, human.Symbol%in% sb_homo$SYMBOL)  
sb_zmh = unique(c(sb_zm$zfin_id_symbol,sb_zh$zfin_id_symbol))


# #social behaviour -------------------------------------------------------

for(i in 1:8){
  deg=get(paste0("deg.cpm",i))
  deg = deg %>% filter(zfin_id_symbol %in% sb_zmh) %>% filter(id %in% b.ft.id) %>% 
    filter(zfin_id_symbol %in% prot_id_map$SYMBOL)
  name = paste0("sb.deg.cpm",i)
  assign(name,deg)
}


#6 dpf
sb.deg.5 = sb.deg.cpm5 %>% filter(zfin_id_symbol %in% sb_zmh) %>% filter(FDR < 0.05)%>% filter(abs(logFC) > log2(1.5))
#13 dpf
sb.deg.6 = sb.deg.cpm6 %>% filter(zfin_id_symbol %in% sb_zmh) %>% filter(FDR < 0.05)%>% filter(abs(logFC) > log2(1.5))
#120 dpf
sb.deg.120 = sb.deg.cpm7 %>% filter(zfin_id_symbol %in% sb_zmh) %>% filter(FDR < 0.05)%>% filter(abs(logFC) > log2(1.5))
#120 dpf after LD
sb.deg.LD = sb.deg.cpm8 %>% filter(zfin_id_symbol %in% sb_zmh) %>% filter(FDR < 0.05)%>% filter(abs(logFC) > log2(1.5))

#DEGs in 6 or 13 or 120
sb.613120= unique(c(sb.deg.5$zfin_id_symbol,sb.deg.6$zfin_id_symbol,sb.deg.120$zfin_id_symbol))
#DEGs in 6 or 13 or 120 or LD
sb.613120ld= unique(c(sb.deg.5$zfin_id_symbol,sb.deg.6$zfin_id_symbol,sb.deg.120$zfin_id_symbol,sb.deg.LD$zfin_id_symbol))

#table for human homologs
sig.sb_zh = filter(zh.all, zfin_id_symbol %in% sb.613120)
sig.sb_zh= na.omit(sig.sb_zh)
#write.csv(sig.sb_zh,"./data/sig.p3.sb_zh.csv")


#heatmap and venn diagrams

for(i in 1:8){
  deg=get(paste0("sb.deg.cpm",i))
  deg = deg %>% filter(zfin_id_symbol %in% sb.613120) 
  name = paste0("d6tom4.sb.deg.cpm",i)
  assign(name,deg)
}


###d6, 13, 120 dpf wt and star.bPAC+
ht.mt = data.frame(d6tom4.sb.deg.cpm5$logFC,
                   d6tom4.sb.deg.cpm6$logFC,
                   d6tom4.sb.deg.cpm7$logFC)
row.names(ht.mt)=d6tom4.sb.deg.cpm5$zfin_id_symbol

ht.FDR= data.frame(d6tom4.sb.deg.cpm5$FDR,
                   d6tom4.sb.deg.cpm6$FDR,
                   d6tom4.sb.deg.cpm7$FDR)
row.names(ht.FDR)=d6tom4.sb.deg.cpm5$zfin_id_symbol

ht.FDR[ht.FDR <= 0.05] <-"*"
ht.FDR[ht.FDR > 0.05] <-NA
ht.FDR[abs(ht.mt) < log2(1.5)] <-NA



# #heat map ---------------------------------------------------------------


library (gplots)
h.breaks=seq(-3,3,0.01)
#mycol <- colorpanel(n=length(breaksList)-1,low='blue', mid = 'white', high='firebrick3')
#mycol <- colorpanel(1000,"blue","white","red")
mycol = colorRampPalette(c("#313695","#4575B4","#74ADD1","#ABD9E9" ,"#E0F3F8","#FFFFBF" ,"#FEE090" ,"#FDAE61", 
                           "#F46D43","#D73027","#A50026" ))(n=600)

st_sig_tx_1= as.matrix(ht.mt)

tiff("./figures/s.Fig4a_sb.6to120.4.5x22.tiff", height = 22, width = 4.5, res = 1200, compression = "lzw", units = "cm")
try(heatmap.2((st_sig_tx_1), col=mycol, scale='none',
               dendrogram = "row",Colv = F,breaks=h.breaks,
               #main=c("social behaviour, log2FC"),
               density='density', trace = 'none',
               na.color = F,
               offsetRow = 0,
               offsetCol = c(0,0.2),
               cexCol = 1,
               labCol = c("d6","d13","d120"),
               cellnote=ht.FDR,
               notecex=1,
               notecol="black",
               margins = c(5, 5),
               #key.title = "log2FC",
               keysize = 1, symkey=F),silent = T)

dev.off()




#Gene expression profile of 
# #oxt regulatory elements ------------------------------------------------

for(i in 1:8){
  deg=get(paste0("deg.cpm",i))
  deg = deg %>% filter(zfin_id_symbol %in% oxt_reg) %>% filter(id %in% b.ft.id) %>% 
    filter(zfin_id_symbol %in% prot_id_map$SYMBOL)
  name = paste0("oxtr.deg.cpm",i)
  assign(name,deg)
}

#6 dpf
oxtr.deg.5 = oxtr.deg.cpm5 %>% filter(zfin_id_symbol %in% oxt_reg) %>% filter(FDR < 0.05)
#13 dpf
oxtr.deg.6 = oxtr.deg.cpm6 %>% filter(zfin_id_symbol %in% oxt_reg) %>% filter(FDR < 0.05)
#120 dpf
oxtr.deg.120 = oxtr.deg.cpm7 %>% filter(zfin_id_symbol %in% oxt_reg) %>% filter(FDR < 0.05)
#120 dof after LD
oxtr.deg.LD = oxtr.deg.cpm8 %>% filter(zfin_id_symbol %in% oxt_reg) %>% filter(FDR < 0.05)

#DEGs in 6 or 13 or 120
oxtr.613120= unique(c(oxtr.deg.5$zfin_id_symbol,oxtr.deg.6$zfin_id_symbol,oxtr.deg.120$zfin_id_symbol))
#DEGs in 6 or 13 or 120 or LD
oxtr.613120ld= unique(c(oxtr.deg.5$zfin_id_symbol,oxtr.deg.6$zfin_id_symbol,oxtr.deg.120$zfin_id_symbol,oxtr.deg.LD$zfin_id_symbol))



# ## tables for heatmap(s.Fig4b) ---------------------------------------------------


for(i in 1:8){
  deg=get(paste0("oxtr.deg.cpm",i))
  deg = deg %>% filter(zfin_id_symbol %in% oxtr.613120) 
  name = paste0("d6.oxtr.deg.cpm",i)
  assign(name,deg)
}


###star:bPAC+
ht.mt = data.frame(d6.oxtr.deg.cpm5$logFC,d6.oxtr.deg.cpm6$logFC,d6.oxtr.deg.cpm7$logFC)
row.names(ht.mt)=d6.oxtr.deg.cpm5$zfin_id_symbol

ht.FDR= data.frame(d6.oxtr.deg.cpm5$FDR,d6.oxtr.deg.cpm6$FDR,d6.oxtr.deg.cpm7$FDR)
row.names(ht.FDR)=d6.oxtr.deg.cpm5$zfin_id_symbol

ht.FDR[ht.FDR <= 0.05] <-"*"
ht.FDR[ht.FDR > 0.05] <-NA
ht.FDR[abs(ht.mt) < log2(1.5)] <-NA

library (gplots)
h.breaks=seq(-3,3,0.01)
#mycol <- colorpanel(n=length(breaksList)-1,low='blue', mid = 'white', high='firebrick3')
#mycol <- colorpanel(1000,"blue","white","red")
mycol = colorRampPalette(c("#313695","#4575B4","#74ADD1","#ABD9E9" ,"#E0F3F8","#FFFFBF" ,"#FEE090" ,"#FDAE61", 
                           "#F46D43","#D73027","#A50026" ))(n=600)


st_sig_tx_1= as.matrix(ht.mt)


# #heat map ---------------------------------------------------------------


tiff("./figures/s.Fig4b_heatamp.oxt.regulators.tiff", height = 12, width = 4.5, res = 1200, compression = "lzw", units = "cm")

try(heatmap.2((st_sig_tx_1), col=mycol, scale='none',
               dendrogram = "row",Colv = F,breaks=h.breaks,
               #main=c("oxt_regulation, log2FC"),
               density='density', trace = 'none',
               na.color = F,
               offsetRow = 0,
               offsetCol = c(0,0.2),
               cexCol = 1,
               labCol = c("d6","d13","d120"),
               cellnote=ht.FDR,
               notecex=1,
               notecol="black",
               margins = c(5, 5),
               #key.title = "log2FC",
               keysize = 0.8 ,symkey=F ),silent = T)

dev.off()


# Chea3_TF_enrichment -----------------------------------------------------

#input

zh.sb.deg = sig.sb_zh$human.Symbol

sb.deg.up = unique(c(sb.deg.5[c(which(sb.deg.5$logFC >0)),8],
                     sb.deg.6[c(which(sb.deg.6$logFC >0)),8],
                     sb.deg.120[c(which(sb.deg.120$logFC >0)),8])) 
sb.deg.up=na.omit(sb.deg.up)
zh.sb.deg.up = filter(sig.sb_zh,zfin_id_symbol %in% sb.deg.up)
zh.sb.deg.up = unique(zh.sb.deg.up$human.Symbol)

sb.deg.dn = unique(c(sb.deg.5[c(which(sb.deg.5$logFC <0)),8],
                  sb.deg.6[c(which(sb.deg.6$logFC <0)),8],
                  sb.deg.120[c(which(sb.deg.120$logFC <0)),8]))
sb.deg.dn=na.omit(sb.deg.dn)
zh.sb.deg.dn = filter(sig.sb_zh,zfin_id_symbol %in% sb.deg.dn)
zh.sb.deg.dn = unique(zh.sb.deg.dn$human.Symbol)

#write(zh.sb.deg.up,"./outputs/sb.up.human_homologs.txt",sep = "\n")
#write(zh.sb.deg.dn,"./outputs/sb.dn.human_homologs.txt",sep = "\n")

###chea3_results
#Chea3 API: https://maayanlab.cloud/chea3/
# run ChIP-X Enrichment analysis version 3

library(httr)
library(jsonlite)

url = "https://maayanlab.cloud/chea3/api/enrich/"
encode = "json"
payload = list(query_name = "sb_gene", gene_set =zh.sb.deg )

#POST to ChEA3 server
response = POST(url = url, body = payload, encode = encode)
json = content(response, "text")  

chea3_sb = fromJSON(json)
overlapped.number= str_count(chea3_sb$`Integrated--meanRank`$Overlapping_Genes,",")+1
chea3_sb.int.mean.rank = cbind(chea3_sb$`Integrated--meanRank`,overlapped.number)

write.csv(chea3_sb.int.mean.rank,"./outputs/s.fig.4E.chea3_sb_results_Integrated--meanRank.csv")  
  

chea3 = data.frame(chea3_sb$`Integrated--meanRank`)
chea3.sb50= chea3 %>% filter(as.numeric(Rank) <= 50) 
chea3.dn.deg = chea3 %>% filter(as.numeric(Rank) <= 50) %>%
  filter(TF %in% c(zh.all$human.Symbol[c(which(zh.all$zfin_id_symbol %in% deg_all.dn$zfin_id_symbol))]))

library(forcats) 

top.tf.sb=chea3.dn.deg$TF
top.tf.sb.rank= paste0("rank:",chea3.dn.deg$Rank)
top.tf.sb.score= log2(1632/as.numeric(chea3.dn.deg$Score))

topTF.zh_sb=filter(ft.deg.cpm.dn5, zfin_id_symbol %in% unique(c(filter(zh.all,human.Symbol %in% top.tf.sb)$zfin_id_symbol)))$zfin_id_symbol
topTF.hz_sb=unique(c(filter(zh.all,zfin_id_symbol %in% topTF.zh_sb)$human.Symbol))


tf.bp.data = data.frame("TF"=top.tf.sb[1:10],"rank"=top.tf.sb.rank[1:10],"score"=top.tf.sb.score[1:10])

col= c()
for (i in 1:10) {
  if (tf.bp.data$TF[i] %in% topTF.hz_sb) {
    col= c(col,"forestgreen")
  }else {
    col= c(col,"grey")
  }
}
tf.bp.data[,"col"]=col
tf.bp.data=mutate(tf.bp.data,TF = fct_reorder(TF, score))

##barplot

sb.tf.p= ggplot(tf.bp.data, aes(x=TF, y=score))+
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


tiff("./figures/s.Fig4E_sb-.tf.tiff", width = 12, height = 12,units = "cm", res = 1200, compression = "lzw")
plot(sb.tf.p)
dev.off()

