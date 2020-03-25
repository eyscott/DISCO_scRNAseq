# DISCO: scRNAseq
DISCO-Digital Microfluidic Isolation of Single Cells for - Omics
The laser lysis capture method implemented in DISCO offers spatially and temporally resolved capture of single cells, while Digital Microfluidics (DMF) is capable of retrieving this lysate and permitting flexible downstream analyses, such as scRNA-seq.
Each cell has a barcode incorporated into Read 1, making the cell's transcriptome accountable and traceable through the pipeline. This barcode includes a unique cell barcode (12bp) and a UMI (8bp), the first part of our bioinformatic pipeline deals with parsing these.

## General Worlflow
The provided scripts go through mapping, parsing fasta headers (barcodes and UMIs), generating count data, making figures and performing differential gene expression (DGE):
1. Append Read1 barcode information to Read2 fasta header 
2. Align with Star Aligner 
3. Collapse redundant UMIs 
4. Parse and separate unique cell barcodes 
5. Gather gene count data (FeatureCounts)

RFigureScripts: TPM normalization, figure generation, DGE with EdgeR, UMAP implementation 

DependentScripts: FeatureCounts and barcode demuliplexing scripts.

## Publication
pending:  
*Lamanna, J.*,*Scott, E.Y.*, *Edwards, H.*, Chamberlain, M.D., Dryden, M.D.M., Peng, J., Mair, B., Lee, A., Sklavounos, A.A., Abbas, F., Moffat, J. & A.R. Wheeler. "Digital Microfluidic Isolation of Single Cells for - Omics" 2020. 
  *Lamanna, J.,Scott, E.Y., Edwards, H. are Co-first authors*

### Authors of scripts
Erica Y. Scott: Scripts 1-5 and Figure scripts 1-7  
Harrison Edwards: python barcode parser scripts within DependentScripts/ 
