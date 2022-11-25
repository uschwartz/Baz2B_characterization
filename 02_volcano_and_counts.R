setwd("~/Analysis/Analysis_sync/S011-BRG1/publish")
load("/Users/admin///Library/Mobile Documents/com~apple~CloudDocs/Organisation_NAC/Cooperate_Design/colors/nac_palette_extend.rda")
palette(nac_palette_extend)

res.table<-read.delim("/Volumes/PromisePegasus/Projects/14_Baz2b_Microarray/analysis/01_conventional/results/res.baz2b_del.txt")

head(res.table)

#upreagulated types 0.05
up.05.types<-table(subset(res.table, adj.P.Val<0.05 & logFC > 0)$biotype)
#downreagulated types 0.05
down.05.types<-table(subset(res.table, adj.P.Val<0.05 & logFC < 0)$biotype)


#upreagulated types 0.1
up.1.types<-table(subset(res.table, adj.P.Val<0.1 & logFC > 0)$biotype)
#downreagulated types 0.1
down.1.types<-table(subset(res.table, adj.P.Val<0.1 & logFC < 0)$biotype)


#sort
names<-unique(c(names(up.1.types),names(down.1.types)))

type.table<-rbind(up.05.types[names], down.05.types[names], 
      up.1.types[names], down.1.types[names])

colnames(type.table)<-names

type.table[is.na(type.table)]<-0

rownames(type.table)<-c("up_0.05","down.0.05","up_0.1", "down.0.1")
write.csv(type.table, "regulated_transcripts.csv")


## volcano plot

library(ggplot)
library(ggrepel)

#### 0.1 cutoff

res.table$signif_0.1<-"no change"
res.table$signif_0.1[with(res.table,(adj.P.Val<0.1 & logFC > log2(1.5)))]<-"up"
res.table$signif_0.1[with(res.table,(adj.P.Val<0.1 & logFC < -log2(1.5)))]<-"down"
res.table$signif_0.1<-factor(res.table$signif_0.1, 
                             levels = c("down","no change","up"))

# plot adding up all layers we have seen so far
g<-ggplot(data=res.table, aes(x=logFC, y=-log10(adj.P.Val), col=signif_0.1, label=SYMBOL)) +
    geom_point() + 
    theme_minimal() +
    scale_color_manual(values=c("#357EBD", "#8E8E8D", "#D43F39"))+
    xlab("log2 Fold Change")

#get top 20 p-adjusted value
df.top<-res.table[order(res.table$adj.P.Val)[1:30],]


pdf("volcano_0.1.pdf", width = 7, height = 3.5)
    print(g+geom_text_repel(data = df.top,
                        aes(label = SYMBOL),
                        max.overlaps = 30,
                        size = 2.5,
                        box.padding = unit(0.3, "lines"),
                        point.padding = unit(0.25, "lines")))
dev.off()



###################### 0.05 cutoff ###################
res.table$signif_0.05<-"no change"
res.table$signif_0.05[with(res.table,(adj.P.Val<0.05 & logFC > log2(1.5)))]<-"up"
res.table$signif_0.05[with(res.table,(adj.P.Val<0.05 & logFC < -log2(1.5)))]<-"down"
res.table$signif_0.05<-factor(res.table$signif_0.05, 
                             levels = c("down","no change","up"))

# plot adding up all layers we have seen so far
g<-ggplot(data=res.table, aes(x=logFC, y=-log10(adj.P.Val), col=signif_0.05, label=SYMBOL)) +
    geom_point() + 
    theme_minimal() +
    scale_color_manual(values=c("#357EBD", "#8E8E8D", "#D43F39"))+
    xlab("log2 Fold Change")

#get top 20 p-adjusted value
df.top<-res.table[order(res.table$adj.P.Val)[1:30],]


pdf("volcano_0.05.pdf", width = 7, height = 3.5)
    print(g+geom_text_repel(data = df.top,
                        aes(label = SYMBOL),
                        max.overlaps = 30,
                        size = 2.5,
                        box.padding = unit(0.3, "lines"),
                        point.padding = unit(0.25, "lines")))
dev.off()



