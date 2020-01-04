# TS1015-DESeq2

These data and scripts are supplement to:
Say, T. E., & Degnan, S. M. (2019). Molecular and behavioural evidence that interdependent photo - and chemosensory systems regulate larval settlement in a marine sponge. Molecular Ecology, 00, 1-15. doi:10.1111/mec.15318


#This folder contains the files and scripts used to perform DESeq2 (Love et al., 2014) on the CEL-Seq2 data. 


Script:
"file name: R_Script_DESeq2_TS1015.sh" - R commands used to perform DESeq2 on the CEL-Seq2 dataset




#This script is divided into different parts, outlined below:

Part a) 	import data into R, set up experimental design formula 

Part b)	transform count data (rlog stabilization and variance stabiliazation)

Part c) 	LRT analysis in Deseq2

Part d)	  PCA

Part e) 	heatmap of differentially expressed genes #previously denoted as f

Part f) 	pariwise comparisons between 1 and 9 hpe for each light regime #prev h

Part g) 	plot counts (log2) for genes of interest


INPUT_Files:
All transcriptome data and the counts data (renamed) from this study are available in the NCBI GEO database under accession number GSE130274. The raw counts file (prior to updating sample names) is named as: TS1015_L001_L002_counting_report_tidy_NoERCC_NolncRNAs_NoBact_2016.02.08.txt. 


OUTPUT_Files:
"TS1015_DESeq2_rlogMat_Blind_2016.11.21" - rlog stabilised transformed counts for all genes in this CEL-Seq dataset.

"TS1015_DESeq2_vsdMat_Blind_2016.11.21" - Variance stabilising transformed (vst) counts for all genes expressed in this CEL-Seq dataset.  

The heatmap is in Figure 1d

Pairwise comparisons between 1 and 9 hpe are represented in the venn diagrams in Figure 2e,f

The PCA is in Figure S2.

The plot counts is in Figure 4b.

Reference

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550. doi:10.1186/s13059-014-0550-8
