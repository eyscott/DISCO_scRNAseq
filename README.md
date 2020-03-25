# DISCO_scRNAseq
DISCO-Digital Microfluidic Isolation of Single Cells for - Omics
The laser lysis capture method implemented in DISCO offers spatially and temporally resolved capture of single cells, while Digital Microfluidics (DMF) is capable of retrieving this lysate and permitting flexible downstream analyses, such as scRNA-seq.
Each cell has a barcode incorporated into Read 1, making the cell's transcriptome accountable and traceable through the pipeline. This barcode includes a unique cell barcode (12bp) and a UMI (8bp), the first part of our bioinformatic pipeline deals with parsing these.

## General Worlflow
The provided scripts go through mapping, parsing fasta headers (barcodes and UMIs), generating count data, making figures and performing differential gene expression (DGE):
1. 
2.
3.
4.
5.
FigureScripts: 
DependentScripts: 

## Publication
pending

### Authors of scripts
Erica Y. Scott: Scripts 1-5 and Figure scripts 1-6 
Harrison Edwards: python barcode parser scripts within DependentScripts/ 
