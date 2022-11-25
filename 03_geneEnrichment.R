setwd("~/Analysis/Analysis_sync/S011-BRG1/publish")
load("/Users/admin///Library/Mobile Documents/com~apple~CloudDocs/Organisation_NAC/Cooperate_Design/colors/nac_palette_extend.rda")
palette(nac_palette_extend)

res.table<-read.delim("/Volumes/PromisePegasus/Projects/14_Baz2b_Microarray/analysis/01_conventional/results/res.baz2b_del.txt")

head(res.table)

#upreagulated types 0.05
res.up.05<-subset(res.table, adj.P.Val<0.05 & logFC > 0)
#downreagulated types 0.05
res.down.05<-subset(res.table, adj.P.Val<0.05 & logFC < 0)


df<-data.frame(name=c("up-regulated", "down-regulated"),
           genes=
           c(paste(res.up.05$ENSEMBL, collapse = ","),
           paste(res.down.05$ENSEMBL, collapse = ",")))


write.table(df, "down_and_up_metascape.txt", col.names = F,
            row.names = F, sep="\t", quote = F)


### universe

write.table(res.table$ENSEMBL, "universe.txt",
            col.names = F,
            row.names = F, sep="\n", quote = F)


write.table(res.up.05$ENSEMBL, "up-regulated.txt",
            col.names = F,
            row.names = F, sep="\n", quote = F)

write.table(res.down.05$ENSEMBL, "down-regulated.txt",
            col.names = F,
            row.names = F, sep="\n", quote = F)
