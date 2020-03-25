tbl_mod <- read.table("Gmouse_R_Norm_table.txt", header = T, stringsAsFactors = F)
row.names(tbl_mod) <- tbl_mod[ ,1]
Gmouse_R_N <- tbl_mod[ ,-c(1,28)]
#convert the matrix data from character to integer values, need to for umap input
Gmouse_R_N_t <- t(Gmouse_R_N[ ,-c(16)])
Gmouse_R_NN_t<- as.matrix(type.convert(Gmouse_R_N_t,na.strings = "NA", as.is = FALSE, dec = "."))

library(umap,lib.loc="/home/eyscott/R/x86_64-pc-linux-gnu-library/3.5")
Gmouse_R_umap_t <- umap(Gmouse_R_NN_t,n_epochs=500,n_neighbors=15)
Gmouse_R_umap_t_lay <- as.data.frame(Gmouse_R_umap_t$layout) 
colnames(Gmouse_R_umap_t_lay)<- c('Dim1','Dim2')
Gmouse_R_umap_t_lay$sample <- rownames(Gmouse_R_umap_t_lay)
write.table(Gmouse_R_umap_t_lay, "Gmouse_R_umap.txt",sep="\t",row.names = F)
