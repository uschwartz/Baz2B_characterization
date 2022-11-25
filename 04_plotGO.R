setwd("~/Analysis/Analysis_sync/S011-BRG1/publish")
GO.final<-read.csv("Metascape/Enrichment_GO/_FINAL_GO.csv")


GOred<-GO.final[,c("X_LogP_down.regulated","X_LogP_up.regulated","GO","Description",
            "GROUP_ID", "X.GeneInGOAndHitList")]

split(GOred, GOred$GROUP_ID)


GOdf<-GOred[which(GOred$GO %in% c("GO:0050678","GO:0098609","GO:0030198",
  "GO:0048880","GO:0000904","GO:0050865","GO:0001934")),]

library(ggplot2)

down.df<-GOdf[,-2]
colnames(down.df)[1]<-"log.p.value"
down.df$regulation<-"down-regulated"

up.df<-GOdf[,-1]
colnames(up.df)[1]<-"log.p.value"
up.df$regulation<-"up-regulated"

final.df<-rbind(up.df,down.df)

final.df$Description<-factor(final.df$Description,levels =rev( 
                    unique(final.df[order(final.df$log.p.value),"Description"])))


g<-ggplot(final.df,aes(x=regulation, y=Description, fill=c(-1)*log.p.value,
                    size=c(-1)*log.p.value))+
    geom_point(colour="black",pch=21)+theme_bw()+
    scale_fill_gradientn(colours=c("white",RColorBrewer::brewer.pal(4,"Oranges")))+
    scale_size(range = c(0.1,5))

pdf("Gene_enrichment_analysis.pdf", height = 2, width = 5)
    print(g)
dev.off()

