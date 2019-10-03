# TS1015-DESeq2

#This folder contains the files and scripts used to perform DESeq2 on the CEL-Seq2 data.

Script:
"file name: R_Script_DESeq2_TS1015.sh" - R commands used to perform DESeq2 on the CEL-Seq2 dataset

This script is divided into different parts, outlined below:

Part a) 	import data into R, set up experimental design formula 

Part b)	transform count data (rlog stabilization and variance stabiliazation)

Part c) 	LRT analysis in Deseq2

Part d)	  PCA

Part e) 	heatmap of differentially expressed genes #previously denoted as f

Part f) 	pariwise comparisons between 1 and 9 hpe for each light regime #prev h

Part g) 	plot counts (log2) for genes of interest

OUTPUT_Files:
"TS1015_DESeq2_rlogMat_Blind_2016.11.21" - rlog stabilised transformed counts for all genes in this CEL-Seq dataset.

"TS1015_DESeq2_vsdMat_Blind_2016.11.21" - Variance stabilising transformed (vst) counts for all genes expressed in this CEL-Seq dataset.  

The heatmap is in Figure 1d

Pairwise comparisons between 1 and 9 hpe are represented in the venn diagrams in Figure 2e,f

The PCA is in Figure S2.

The plot counts is in Figure 4b.
