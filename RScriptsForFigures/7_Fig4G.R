QLF_heatmap<- merge(qlf_heatmap_mat,tagsbysample_qlf, by = "Gene_id")
QLF_heatmap_1 <- QLF_heatmap[ ,c(7:19)]
rownames(QLF_heatmap_1)<- QLF_heatmap_1$Gene.name
QLF_heatmap_2<- QLF_heatmap_1[ ,-c(1)]
QLF_heatmap_2mat <- as.matrix(type.convert(QLF_heatmap_2,na.strings = "NA", as.is = FALSE, dec = "."))

#cluster the rows
hr <- hclust(as.dist(1-cor(t(QLF_heatmap_2mat), method="pearson")),
             method="average")
hc <- hclust(as.dist(1-cor(QLF_heatmap_2mat, method="spearman")), method="average")


lmat = rbind(c(0,4),c(0,3),c(2,1))
lwid = c(0.5,4)
lhei = c(1,0.5,4)

library(gplots)
library(viridis)
library(ggridges)

write.table(QLF_heatmap_2mat,"DGE_QLF_heatmap_input.txt")
QLF_heatmap_2 <- read.table("DGE_QLF_heatmap_input.txt")
QLF_heatmap_2mat <- as.matrix(type.convert(QLF_heatmap_2,na.strings = "NA", as.is = FALSE, dec = "."))

png(filename = "Red_QLF_heatmap.png",
    width = 480, height = 800, units = "px", pointsize = 20,
    bg = "white")
par(mar=c(9,6,6,3)+0.2, pin=c(0,0)) 
col_breaks <-c(seq(0,2.4,0.2),seq(3,5,0.5),6,7.5,10,12.5,15,17.5,20,50,75,100,150,200,300,400,500)
heatmap.2(QLF_heatmap_2mat,    
          trace="none", 
          labRow = FALSE,
          margins =c(8,8),  
          col=viridis(32),       
          breaks=col_breaks, 
          labCol = c("Single B16-Pure 1","Single B16-Pure 2","Single B16-Pure 3","Single B16-Pure 4",
                     "Single B16-Pure 5","Single B16-Pure 6","Single B16-Pure 7","Single B16-Pure 8",
                     "Single B16-Mix 1","Single B16-Mix 2","Single B16-Mix 3","Single B16-Mix 4"),   
          cexCol = 1,
          dendrogram="none",   
          Colv=F,
          offsetCol=0,
          adjCol=c(1,1.01),
          adjRow=c(0.3,0),
          Rowv=as.dendrogram(hr),
          hclustfun = hclust,
          keysize = 1,
          key.xlab = "CPM",
          key.ylab = "density",
          key.title = NA,
          key.par=list(mar=c(0.5,1.3,2,6.7)), #BLTR
          densadj = 0.5,
          density.info="density",
          denscol = "white",
          lmat = lmat,
          lwid = lwid,
          lhei = lhei)                
dev.off()  

hmap_order <- data.frame(QLF_heatmap_2mat[rev(hr$labels[hr$order]), hc$labels[hc$order]])
write.csv(hmap_order ,"QLF_heatmap.csv")
