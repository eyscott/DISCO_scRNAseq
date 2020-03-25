#stuck the ENSEMBL gene id versions into Biomart to get further metadata
biomart_R <- read.csv("R_biomart_output.txt",header=T)

geneR_exp_NORM_meta <- merge(geneR_exp_NORM,biomart_R,by.x="Gene_id",by.y="Gene.stable.ID.version")
geneR_exp_NORM_meta_unique <- distinct(geneR_exp_NORM_meta, Gene_id, .keep_all= TRUE)
write.table(geneR_exp_NORM_meta_unique, "Red_NORM_meta.txt")

##mold for figures
geneR_exp_NORM_meta_unique_red <- geneR_exp_NORM_meta_unique[ ,c(1:16,21,22)]
geneR_exp_NORM_meta_unique_melt <- gather(geneR_exp_NORM_meta_unique_red,"variable", "value",2:16)

library(RColorBrewer)
n <-40
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

##Making gene content pie####
library(ggplot2)
pdf(file='Red_gene_pie.pdf', width=10, height=6,bg="white")
ggplot(geneR_exp_NORM_meta_unique_melt,aes(x=factor(1),weight=value,fill=Gene.type)) + 
  geom_bar(width=1,alpha = 1) + xlab(" ") + 
  ylab(" ") + coord_polar("y") + 
  guides(fill=guide_legend(title="RNA type")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 10)) +
  theme(axis.text = element_text(colour="black", size = 8)) +
  scale_fill_manual(values = col_vector) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())
dev.off()
