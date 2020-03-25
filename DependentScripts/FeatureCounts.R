library(Rsubread)
files_G<- c("G1_S4_L001_human_starAligned.sortedByCoord.out_unique_only1.sam","G1_S4_L001_human_starAligned.sortedByCoord.out_unique_only2.sam","G1_S4_L001_human_starAligned.sortedByCoord.out_unique_only3.sam","G1_S4_L001_human_starAligned.sortedByCoord.out_unique_only4.sam","G1_S4_L001_human_starAligned.sortedByCoord.out_unique_only5.sam","Gmix_S7_L001_human_starAligned.sortedByCoord.out_unique_only1.sam","Gmix_S7_L001_human_starAligned.sortedByCoord.out_unique_only2.sam","Gmix_S7_L001_human_starAligned.sortedByCoord.out_unique_only3.sam","Gmix_S7_L001_human_starAligned.sortedByCoord.out_unique_only4.sam","Gmix_S7_L001_human_starAligned.sortedByCoord.out_unique_only5.sam")
FC_G <- featureCounts(files = files_G, annot.ext = '/home/eyscott/project/eyscott/eyscott/RG_scRNA/human_genome/gencode.v26.primary_assembly.annotation.gtf', isGTFAnnotationFile = TRUE)
write.table(FC_G$stat,"FC_G_stat.txt", quote=F, sep = "\t")
write.table(FC_G$counts,"FC_G_counts.txt", quote=F, sep = "\t")
write.table(FC_G$annotation,"FC_G_annot.txt", quote=F, sep = "\t")

files_R <- files_R
FC_R <- featureCounts(files = files_R, annot.ext = '/home/eyscott/project/eyscott/eyscott/RG_scRNA/mouse_genome/gencode.vM15.primary_assembly.annotation.gtf', isGTFAnnotationFile = TRUE)
write.table(FC_R$stat,"FC_R_stat.txt", quote=F, sep = "\t")
write.table(FC_R$counts,"FC_R_counts.txt", quote=F, sep = "\t")
write.table(FC_R$annotation,"FC_R_annot.txt", quote=F, sep = "\t")
