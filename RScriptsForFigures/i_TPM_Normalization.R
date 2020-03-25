###TPM Normalization###
#Get gene length
geneR_L<- read.table("FC_R_annot.txt", header = T, stringsAsFactors = F)
#get number of reads
R_stats <- as.data.frame(t(read.table("FC_R_stat.txt", header = T, stringsAsFactors = F,row.names = NULL)))
R_stats_red <- as.data.frame(R_stats[-c(1:2) ,1])
colnames(R_stats_red) <- "Assigned"
row.names(R_stats_red) <- c("R1_S","R2_S","R3_S","R4_S","R5_sPos","R1_S2","R2_S2","R3_S2","R4_S2","R5_sPos2","R1_mix","R2_mix","R3_mix","R4_mix","R5_mPos")
row.names(R_stats_red) -> R_stats_red$variable

###Get expression data
library(tidyr)
library(plyr)
library(dplyr)
geneR_exp <- read.table("FC_R_counts.txt", header = T, stringsAsFactors = F)
colnames(geneR_exp) <- c("R1_S","R2_S","R3_S","R4_S","R5_sPos","R1_S2","R2_S2","R3_S2","R4_S2","R5_sPos2","R1_mix","R2_mix","R3_mix","R4_mix","R5_mPos")
row.names(geneR_exp) -> geneR_exp$Gene_id
geneR_exp_melt <- gather(geneR_exp,"variable", "value", 1:15)
##normalize for gene length
geneR_exp_melt_L <- merge(geneR_exp_melt,geneR_L, by.x="Gene_id", by.y="GeneID")
geneR_exp_melt_L$Length_kb <- (geneR_exp_melt_L$Length)/1000
geneR_exp_melt_L$value_L_norm <- (geneR_exp_melt_L$value/geneR_exp_melt_L$Length_kb)
#now normalize for library depth
geneR_exp_melt_L_presum<- ddply(geneR_exp_melt_L, c("variable"), summarise,
                                Count_sum = sum(value_L_norm))
geneR_exp_melt_L_presum$RPM <- ((geneR_exp_melt_L_presum$Count_sum)/1000000)
geneR_exp_melt_L_RPM <- merge(geneR_exp_melt_L, geneR_exp_melt_L_presum, by="variable")
geneR_exp_melt_L_RPM$value_L_RPM_norm <- (geneR_exp_melt_L_RPM$value_L_norm/geneR_exp_melt_L_RPM$RPM)

##get rid of some of the temp columns before making gene matrix
geneR_exp_melt_L_RPM_red <- geneR_exp_melt_L_RPM[ ,c("Gene_id", "variable", "value_L_RPM_norm")]
##make the gene matrix again with new normalized values
geneR_exp_NORM <- spread(geneR_exp_melt_L_RPM_red,variable,value_L_RPM_norm)
rownames(geneR_exp_NORM) <- geneR_exp_NORM$Gene_id
geneR_exp_NORM_mat <- as.matrix(geneR_exp_NORM[ ,-1])

setwd('/Users/erica/Desktop/scRNA_analysis_mat/RG_analysis_reps')
write.table(geneR_exp_NORM_mat, "Red_NORM_matrix.txt")
