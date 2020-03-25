Green_mat <- read.table("Green_NORM_matrix.txt", header=T)
Red_mat <- read.table("Red_NORM_matrix.txt", header=T)
#attach gene names to both
rownames(Green_mat)->Green_mat$Gene_id
rownames(Red_mat)->Red_mat$Gene_id
#import gene name
biomart_G <- read.table("G_biomart_output.txt",header=T)
biomart_R <- read.csv("R_biomart_output.txt",header=T)
G_names <- biomart_G[ ,c(1,6)]
R_names <- biomart_R[ ,c(1,6)]
##Now mapped Green human cells to mouse genome to get the same units for genome coordinates
Gmouse_exp <- read.table("Gmouse_NORM_matrix.txt", header=T)
rownames(Gmouse_exp)->Gmouse_exp$Gene_id
biomart_Gmouse <- read.csv("Gmouse_biomart_output.txt",header=T)
Gmouse_names <- biomart_Gmouse[ ,c(1,6)]
Gmouse_mat_names1<- merge(Gmouse_exp,Gmouse_names,by.x='Gene_id',by.y="GenestableIDversion")
##merge red and green matrices
Gmouse_and_R <- merge(Red_mat_names1,Gmouse_mat_names1,by='Gene_id') 
Gmouse_and_R_unique <- distinct(Gmouse_and_R, Gene_id, .keep_all= TRUE) #48315 transcripts
write.table(Gmouse_and_R_unique,"Gmouse_R_Norm_table.txt")

##umap analysis####
Gmouse_R_umap <-read.table("Gmouse_R_umap.txt", header=T, stringsAsFactors = F)
colnames(Gmouse_R_umap)<-c("Dim1","Dim2","Sample")

Gmouse_R_umap$cond[grepl( "R[1-5]" , Gmouse_R_umap$Sample)]<-"B16"
Gmouse_R_umap$cond[grepl( "R[1-5]_mix" , Gmouse_R_umap$Sample)]<-"B16mix"
Gmouse_R_umap$cond[grepl( "R5_sPos" , Gmouse_R_umap$Sample)]<-"B16pos"
Gmouse_R_umap$cond[grepl( "R5_mPos" , Gmouse_R_umap$Sample)]<-"posMix"
Gmouse_R_umap$cond[grepl( "G[1-5]" , Gmouse_R_umap$Sample)]<-"U87"
Gmouse_R_umap$cond[grepl( "G[1-5]mouse_mix" , Gmouse_R_umap$Sample)]<-"U87mix"
Gmouse_R_umap$cond[grepl( "G5mouse_sPos" , Gmouse_R_umap$Sample)]<-"U87pos"
Gmouse_R_umap$cond[grepl( "G5mouse_mPos" , Gmouse_R_umap$Sample)]<-"posMix"
Gmouse_R_umap_sample.names<-Gmouse_R_umap$Sample

pdf(file="Gmouse_red_Umap.pdf",width=7, height=4)
ggplot(Gmouse_R_umap, aes(Dim1,Dim2,colour=cond)) +
  geom_point(size = 5) +
  scale_colour_manual("Cell Type", values = c("#f77d72","#f01d0a","#871f18","#cfc404","#42f73b","#1ec217","#086b04"),
                      labels=c("Single B16-Pure","Single B16-Mix","Bulk B16-Pure","Bulk Mix","Single U87-Pure","Single U87-Mix","Bulk U87")) +
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=16, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(legend.title.align = 0.5) +
  theme(axis.text.y = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(colour="black", size = 12, hjust=1)) +
  theme(axis.title = element_text(colour="black", size = 16)) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, size = 2)) +
  guides(shape = guide_legend(override.aes = list(size=4))) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()
