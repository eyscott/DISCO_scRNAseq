library(edgeR)
setwd('/Users/erica/Desktop/scRNA_analysis_mat/RG_analysis_reps')
#get raw data again
geneR_exp <- read.table("FC_R_counts.txt", header = T, stringsAsFactors = F)
colnames(geneR_exp) <- c("R1_S","R2_S","R3_S","R4_S","R5_sPos","R1_S2","R2_S2","R3_S2","R4_S2","R5_sPos2","R1_mix","R2_mix","R3_mix","R4_mix","R5_mPos")
row.names(geneR_exp) -> geneR_exp$Gene_id
##we are going to attach gene type again and remove junk RNA detected
biomart_R <- read.csv("R_biomart_output.txt",header=T)
geneR_exp_meta <- merge(geneR_exp,biomart_R,by.x="Gene_id",by.y="Gene.stable.ID.version")
geneR_exp_meta_unique <- distinct(geneR_exp_meta, Gene_id, .keep_all= TRUE) #48315
geneR_exp_meta_unique_RNA<-subset(geneR_exp_meta_unique,Gene.type %in% c('lncRNA','miRNA','processed_pseudogene','protein_coding')) #reduces to 37621, reduced by 10694
#rownames(geneR_exp_meta_unique_RNA)<-geneR_exp_meta_unique_RNA$Gene.name
geneR_exp_meta_unique_RNA_singles <- geneR_exp_meta_unique_RNA[ ,c(1:5,7:10,12:15,21)]
rownames(geneR_exp_meta_unique_RNA_singles)<-geneR_exp_meta_unique_RNA_singles$Gene_id
geneR_exp_meta_unique_RNA_singles_mat<-geneR_exp_meta_unique_RNA_singles[ ,-c(1,14)]
geneR_exp_geneNames <- geneR_exp_meta_unique_RNA_singles[ ,c('Gene_id','Gene.name')]
#start setting up edgeR exp
group <- factor(c(1,1,1,1,1,1,1,1,2,2,2,2))
edgeR_counts <- DGEList(counts=geneR_exp_meta_unique_RNA_singles_mat,group=group)
design <- model.matrix(~group)
#calculate normalization factors using TMM
edgeR_counts_norm = calcNormFactors(geneR_exp_meta_unique_RNA_singles_mat, method='TMM')
#calculate dispersion
edgeR_disp<- estimateCommonDisp(geneR_exp_meta_unique_RNA_singles_mat,group=group,lib.size=edgeR_counts_norm, tol=1e-06,
                                rowsum.filter=5, verbose=FALSE)
dim(edgeR_counts)
#[1] 37621    12
#testing using GLM
y <- estimateDisp(edgeR_counts,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
qlf_topTags <- topTags(qlf, n=37621)
#get individualcpm normalized values for each sample to make a heatmap
qlf_topTags_ind <-rownames(topTags(qlf, n=37621)$table)
tagsbysample_qlf<-as.data.frame(cpm(edgeR_counts)[qlf_topTags_ind, order(edgeR_counts$samples$group)])
tagsbysample_qlf$Gene_id <-rownames(tagsbysample_qlf)
#retrieve genes with a significant pvalue
qlf_topTags_id<-subset(qlf_topTags$table, FDR < 0.05) #5300
qlf_topTags_id$Gene_id <-rownames(qlf_topTags_id)
qlf_heatmap_mat <- merge(qlf_topTags_id,geneR_exp_geneNames, by="Gene_id")

#some QC/supplementary plots
plotQLDisp(fit)

pch <- c(0,1)
colors <- rep(c("darkgreen", "red"), 2)
plotMDS(y, col=colors[group], pch=pch[group])
legend("topleft", legend=levels(group), pch=pch, col=colors, ncol=2)

plotMD(y, column=1)
abline(h=0, col="red", lty=2, lwd=2)

plotBCV(y)
