#obtain cumulative TPM for each protein-coding, lncRNA, miRNA and processed-pseudogenes
geneR_exp_NORM_meta_unique_RNA_PROPs_melt_sum<- ddply(geneR_exp_NORM_meta_unique_RNA_stats_melt, c("variable","Gene.type"), summarise,
                                                      Prop_sums = sum(value), Transcript_count=length(Gene.name),
                                                      meanTPM=mean(value))
write.table(geneR_exp_NORM_meta_unique_RNA_PROPs_melt_sum, "Red_RNAproportion_stats.txt")
#Get Gene counts from raw
geneR_exp_melt_noLow <- subset(geneR_exp_melt,value>0.1)
geneR_exp_melt_noLow_counts<- ddply(geneR_exp_melt_noLow, c("variable"), summarise,
                                                      Transcript_count=length(Gene_id))
R_stats_mod <- merge(R_stats_red,geneR_exp_melt_noLow_counts,by='variable')

##scatter of reads vs genes
geneR_exp_NORM_meta_unique_RNA <- subset(geneR_exp_NORM_meta_unique,Gene.type %in% c('lncRNA','protein_coding','processed_pseudogene','miRNA'))
geneR_exp_NORM_meta_unique_RNA_stats <- geneR_exp_NORM_meta_unique_RNA[ ,c(1:16,20:22)]
geneR_exp_NORM_meta_unique_RNA_stats_melt <- gather(geneR_exp_NORM_meta_unique_RNA_stats,'variable','value',2:16)
geneR_exp_NORM_meta_unique_RNA_stats_melt_sum<- ddply(geneR_exp_NORM_meta_unique_RNA_stats_melt, c("variable"), summarise,
                                                      sum = sum(value), gene_L=sum(Transcript.length..including.UTRs.and.CDS.),gene_N=length(Gene.name),
                                                      meanTPM=mean(value))
library(caroline)
geneR_exp_NORM_meta_unique_RNA_PROPs_melt_sum_red <- geneR_exp_NORM_meta_unique_RNA_PROPs_melt_sum[ ,c('variable','Gene.type','Prop_sums')]
R1_S = nv((subset(geneR_exp_NORM_meta_unique_RNA_PROPs_melt_sum_red,variable=='R1_S')$Prop_sums),name=c('lncRNA','miRNA','processed_pseudogene','protein_coding'))
R2_S = nv((subset(geneR_exp_NORM_meta_unique_RNA_PROPs_melt_sum_red,variable=='R2_S')$Prop_sums),name=c('lncRNA','miRNA','processed_pseudogene','protein_coding'))
R3_S = nv((subset(geneR_exp_NORM_meta_unique_RNA_PROPs_melt_sum_red,variable=='R3_S')$Prop_sums),name=c('lncRNA','miRNA','processed_pseudogene','protein_coding'))
R4_S = nv((subset(geneR_exp_NORM_meta_unique_RNA_PROPs_melt_sum_red,variable=='R4_S')$Prop_sums),name=c('lncRNA','miRNA','processed_pseudogene','protein_coding'))
R5_sPos = nv((subset(geneR_exp_NORM_meta_unique_RNA_PROPs_melt_sum_red,variable=='R5_sPos')$Prop_sums),name=c('lncRNA','miRNA','processed_pseudogene','protein_coding'))
R1_S2 = nv((subset(geneR_exp_NORM_meta_unique_RNA_PROPs_melt_sum_red,variable=='R1_S2')$Prop_sums),name=c('lncRNA','miRNA','processed_pseudogene','protein_coding'))
R2_S2 = nv((subset(geneR_exp_NORM_meta_unique_RNA_PROPs_melt_sum_red,variable=='R2_S2')$Prop_sums),name=c('lncRNA','miRNA','processed_pseudogene','protein_coding'))
R3_S2 = nv((subset(geneR_exp_NORM_meta_unique_RNA_PROPs_melt_sum_red,variable=='R3_S2')$Prop_sums),name=c('lncRNA','miRNA','processed_pseudogene','protein_coding'))
R4_S2 = nv((subset(geneR_exp_NORM_meta_unique_RNA_PROPs_melt_sum_red,variable=='R4_S2')$Prop_sums),name=c('lncRNA','miRNA','processed_pseudogene','protein_coding'))
R5_sPos2 = nv((subset(geneR_exp_NORM_meta_unique_RNA_PROPs_melt_sum_red,variable=='R5_sPos2')$Prop_sums),name=c('lncRNA','miRNA','processed_pseudogene','protein_coding'))
R1_mix = nv((subset(geneR_exp_NORM_meta_unique_RNA_PROPs_melt_sum_red,variable=='R1_mix')$Prop_sums),name=c('lncRNA','miRNA','processed_pseudogene','protein_coding'))
R2_mix = nv((subset(geneR_exp_NORM_meta_unique_RNA_PROPs_melt_sum_red,variable=='R2_mix')$Prop_sums),name=c('lncRNA','miRNA','processed_pseudogene','protein_coding'))
R3_mix = nv((subset(geneR_exp_NORM_meta_unique_RNA_PROPs_melt_sum_red,variable=='R3_mix')$Prop_sums),name=c('lncRNA','miRNA','processed_pseudogene','protein_coding'))
R4_mix = nv((subset(geneR_exp_NORM_meta_unique_RNA_PROPs_melt_sum_red,variable=='R4_mix')$Prop_sums),name=c('lncRNA','miRNA','processed_pseudogene','protein_coding'))
R5_mPos = nv((subset(geneR_exp_NORM_meta_unique_RNA_PROPs_melt_sum_red,variable=='R5_mPos')$Prop_sums),name=c('lncRNA','miRNA','processed_pseudogene','protein_coding'))

props <- rbind(R1_S,R2_S,R3_S,R4_S,R5_sPos,R1_S2,R2_S2,R3_S2,R4_S2,R5_sPos2,R1_mix,R2_mix,R3_mix,R4_mix,R5_mPos)
props2 <- cbind(names = rownames(props), props)
props3<-merge(R_stats_mod_num_bind,props2,by.x='variable',by.y='names')
library(scatterpie)
RNA_type<-c('lncRNA','miRNA','processed_pseudogene','protein_coding')
props3$lncRNA <- as.numeric(as.character(props3$lncRNA))
props3$miRNA <- as.numeric(as.character(props3$miRNA))
props3$processed_pseudogene <- as.numeric(as.character(props3$processed_pseudogene))
props3$protein_coding <- as.numeric(as.character(props3$protein_coding))
props3$Assigned <- as.numeric(as.character(props3$Assigned))
props3$Transcript_count <- as.numeric(as.character(props3$Transcript_count))

write.table(props3,"Red_Aligned_geneC_stats.txt")

options("scipen"=100, "digits"=4)
set.seed(80)
props3[, 4:7] <- props3[, 4:7] /100
props3[, 2] <- props3[, 2] /75
props4<- props3[-c(1,4,7,10,13),]
labels <- c('','','','','','','','','*','*')
colours <- c('black','black','black','black','black','black','black','black','blue','blue')
x.labels<-c('1000'=75,'2000'=150,'3000'=225,'4000'=300, '5000'=375, '6000'=450, '7000'=525, '8000'=600,'9000'=675)
x.limits <-seq(1000,9000,1000)
pdf(file='Red_scatterPie.pdf', width=10, height=8,bg="white")
ggplot() +
  geom_scatterpie(aes(x=Assigned, y=Transcript_count, group=variable,r=500),
                  cols=RNA_type, 
                  alpha= 0.8, color=NA, data=props4) +
  coord_fixed() +
  scale_x_discrete(limits = x.limits,labels=x.labels,expand = c(0.1, 0.1)) +
  scale_y_continuous(limits = c(4000,13000),breaks=seq(4000,13000,1000),expand = c(0, 0)) + ylab("Gene count") + xlab("Thousands of Reads") +
  scale_fill_manual(name  ="RNA type",
                    breaks = c('protein_coding','lncRNA','miRNA','processed_pseudogene'),
                    values=my.cols,
                    labels=c("protein-coding","lncRNA","miRNA","processed-pseudogene")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=20, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 18)) +
  theme(axis.title = element_text(colour="black", size = 20, face="bold")) +
  theme(legend.title.align = 0.2) +
  theme(axis.text.y = element_text(colour="black", size = 12,vjust = 1)) +
  theme(axis.text.x = element_text(colour="black", size = 12,angle=35)) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 2)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  geom_text(size=10,aes(x=props4$Assigned, y=props4$Transcript_count,label=labels),hjust=-0.75, vjust=-0.75,show.legend = FALSE)
dev.off()
