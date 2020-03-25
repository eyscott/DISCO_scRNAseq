#consolidate more meta data information
geneR_exp_NORM_meta_unique_chr<- geneR_exp_NORM_meta_unique[ ,c(1:17,21,22)]
geneR_exp_NORM_meta_unique_chr_melt <- gather(geneR_exp_NORM_meta_unique_chr,'variable','value',2:16)

#add chr start to order along chr, then collapse and merge variable column
geneR_exp_NORM_meta_unique_chrFig1<- geneR_exp_NORM_meta_unique[ ,c(1:18,22)]
geneR_exp_NORM_meta_unique_chrFig2 <- geneR_exp_NORM_meta_unique_chrFig1[with(geneR_exp_NORM_meta_unique_chrFig1, order(Chromosome.scaffold.name, Gene.start..bp.)), ]
geneR_exp_NORM_meta_unique_chrFig2_melt <- gather(geneR_exp_NORM_meta_unique_chrFig2,'variable','value',2:16)
ddply_summ_data2 <- ddply(geneR_exp_NORM_meta_unique_chrFig2_melt, .(Chromosome.scaffold.name,Gene_id,Gene.type), summarise, extra = list(as.character(variable,value,Gene.start..bp.)))

chrs_N <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19","X","Y")    
my.cols <- c('#D95F02','#7570B3','#666666','#A6CEE3')

###Making genome coverage plot####
pdf(file='Red_chr_coordnew_2.pdf', width=12, height=6,bg="white")
ggplot(data=subset(ddply_summ_data2,Gene.type %in% c('lncRNA','protein_coding','processed_pseudogene','miRNA') & !is.na(Gene.type)), aes(Chromosome.scaffold.name)) + 
  geom_bar(aes(Chromosome.scaffold.name,group=Gene_id,fill=Gene.type),stat = "count",alpha = 1) + 
  scale_x_discrete(limits = chrs_N, labels = chrs_N) + 
  scale_y_continuous(limits = c(0,4000),breaks=seq(0,4000,1000),expand = c(0, 0)) + ylab("Gene count") + xlab("Chromosome") +
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
  theme(axis.text.x = element_text(colour="black", size = 12)) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 2)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) 
dev.off()
