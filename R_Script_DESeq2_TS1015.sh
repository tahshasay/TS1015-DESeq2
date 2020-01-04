
####################################################################################
# This script is supplement to Say, T. E., & Degnan, S. M. (2019). 
# Interdependent photo- and chemosensory systems regulate larval settlement in a marine sponge. bioRxiv, 519512. doi:10.1101/519512
#
# bioinformatic analysis of Cel-seq2 data using DESeq2 (Love et al., 2014)
#
#Reference
#Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550. doi:10.1186/s13059-014-0550-8
#
# This script was customised based on the DESeq2 manual and vignette (2014 + 2016 v)
# http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
# and from Dave Wheeler's blog at Massey University
# http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
# Bench to bioinformatics + 
#
####################################################################################


getwd() # get working directory
# setwd("~/path/to/working/directory/")
# directory <- "/path/to/counts/directory/"
setwd("/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_R_Scripts_wd_03.00/TS1015_a03.04_DEseq2_DEA_setup_Design-tm.trt.int_res_2016.10.31") ## to change

#check wd. 
getwd()

library("DESeq2") # start R session and load the DESeq2 package
library("ggplot2") # figures

##########################################################################################
# 1. Import raw count data
#
## Name must include Mat (matrix) ie excel file. 
##########################################################################################

TS1015=read.table("/Users/tahshasay/Documents/RAW_COUNT_DATA/TS1015_binary_data_for_R/TS1015_L001_L002_counting_report_tidy_NoERCC_NoBact_2016.02.08.txt", header=TRUE,row.names=1)
#header=TRUE indicates that the first line contains column names and row.names=1 means that the first column should be used as row names.
head(TS1015) # print the first 7 rows of the file
tail(TS1015)
dim(TS1015)  # print the dimensions of the table

# 47411, 43

# --------------------------------------------------------------------------------------
# filtering count data
# Edge R user guide: p 12/ 105. Section "2.6 Filtering"
# DESeq vignette p14. 


# Counts of atleast X in at least Y.
# X = 1 (arbitrary but reasonable, can be changed)
# Y = 7 (size of smallest replicate group)

filteredCounts = TS1015[rowSums(TS1015 >= 1) >= 7,]

head(filteredCounts) # print the first 7 rows of the file
tail(filteredCounts)
dim(filteredCounts)  # print the dimensions of the table

# 28061, 43


#########################################################################################
# 
# 2: Subset the data
#
#########################################################################################


# select relevant samples (for particular analysis)
expTS1015 <- subset(filteredCounts, select=c(1:43)) #### to change
write.csv(as.data.frame(expTS1015),file="expTS1015.csv") #### to change
dim(expTS1015)


#########################################################################################
# 
# 3: Set up the design formula
#
#########################################################################################

# Subsetting the columns of a DESeq data set:
# Create a factor which describes the condition for each sample so that DESeq knows what data each column in the table represents. In the case where multiple libraries have the same condition, they are treated as biological replicates.

Design=data.frame(row.names=colnames(expTS1015),
			treatment=c(rep("drk", 7), rep("lgt",7), rep("nat",7), rep("drk", 7), rep("lgt", 7), rep("nat", 8)), 
			time=c(rep("t1", 21), rep("t9", 22))) # to change

Design


##################
# DESIGN INFO
##################
#Experimental design includes the light regime and time (hours post emergence) for each sample. Light regimes include include natural day-night cycle (nat), constant dark (drk) and constant light (lgt). Time (hours post emergence) include 1 hpe (t1) and 9 hpe (t9). 
#########################################################################################





# create the central data object (structure) in the DESeq package. This is used to store the input values, intermediated calculations and results of an analysis of differential drkression. you can create this object from various types of input: here I use a Matrix of count reads and thus call the function DESeqDataSetFromMatrix construct the data object from the matrix of counts and the metadata table.

# Note: In order to benefit from the default settings of the package, you should put the variable of interest LAST, and make sure the control level is FIRST.
 
 
# To now construct the DESeqDataSet (dds) object from the matrix of counts and the sample information table, we use:
dds<-DESeqDataSetFromMatrix(countData=expTS1015,
					colData=Design,
					design= ~ time + treatment + time:treatment) # to change
	# 	you should put the variable of interest at the end of the formula and make sure the control level is the first level.		 

dds

#print order of factors
colData(dds)$time
dds$time

colData(dds)$treatment
dds$treatment


#########################################################################################
# 
# 4: Set the base level for comparisons (re-level to compare control to treatment)
#
#########################################################################################

# If you call colData - the rows correspond to columns of count data matrix.  The levels in colData are important because they are used in the log calculations.  Reorder treatments so that negative control is the first factor 

# dds$time <- relevel(dds$time, "t0", "t8", "nonind", "sett")
levels(dds$time) 
dds$time <- factor(dds$time, levels=c("t1","t9"))
levels(dds$time) # print the new order of factors

levels(dds$treatment)
dds$treatment <- factor(dds$treatment, levels=c("nat","lgt", "drk"))
levels(dds$treatment) # print the new order of factors


#print order of factors + lists
colData(dds)$time
colData(dds)$treatment

#[1] t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t1 t9 t9 t9 t9 t9 t9 t9 t9 t9 t9 t9 t9 t9
#[35] t9 t9 t9 t9 t9 t9 t9 t9 t9
#Levels: t1 t9


# [1] drk drk drk drk drk drk drk lgt lgt lgt lgt lgt lgt lgt nat nat nat nat nat nat nat drk drk drk drk drk
# [27] drk drk lgt lgt lgt lgt lgt lgt lgt nat nat nat nat nat nat nat nat
# Levels: nat lgt drk

##########################################################################################
###
### gut of the DESeq2 analysis
###
##########################################################################################

dds <- DESeq(dds)
## to change can include test="LRT", reduced= ~1 (to identify all possible DEGs)


# This function will print out a message for the various steps it performs. 
# These are described in more detail in the manual page for DESeq, which can be accessed 
# by typing ?DESeq. Briefly these are: the estimation of size factors (which control for 
# differences in the library size of the sequencing drkeriments), 
# the estimation of dispersion for each gene, and fitting a generalized linear model
# A DESeqDataSet is returned which contains all the fitted information within it, 
# and the following section describes how to extract out results tables of interest from 
# this object.

# look at all of the diff comparisons that were performed during above step
resultsNames(dds)

# output example: "Intercept" "Time_T45_vs_Tbase" "Time_Tend_vs_Tbase" "Genotype_mutant_vs_WT" "TimeT45.Genotypemutant" [6] "TimeTend.Genotypemutant"

# with no arguments to results, the results will be for the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable over the first level. ie. treatment: natural v. light. The text tells you the estimates are of log fold changed logs (trreated/ untreated). 

# res <- results(dds)


########################################################
# print information needed for my excel sheet      # # # 
########################################################

dim(TS1015)
dim(expTS1015)

writeLines(capture.output(dim(TS1015)), "metadata_TS1015_DEA_setup_Design-tm.trt.int_2016.10.31.txt")

writeLines(capture.output(dim(expTS1015)), "metadata_expTS1015_DEA_setup_Design-tm.trt.int_2016.10.31.txt")

writeLines(capture.output(dim(filteredCounts)), "metadata_TS1015_filteredcounts_DEA_setup_Design-tm.trt.int_2016.10.31.txt")


resultsNames(dds)

Design
sessionInfo()

# Capture the screen output into a character vector and use writeLines.
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")


# ---------------------------------------------------
# Save results for use in subsequent parts/ scripts
save(Design, dds, file = "Environment_TS1015_a03.04_DESeq2_DEA_setup_Design-tm.trt.int_2016.10.31.RData")
# --------------------------------------------------- date will change automatically

# load above data in next script:
#__MEs= load(file = "./Environment_TS1015_a03.03_DEA_setup_tm.trt.int_2016.10.26.RData")



# how to save r console output to file automatically
# https://jdegen.wordpress.com/2011/02/09/how-to-save-r-console-output-to-file-automatically/
# say you want to run a series of models and save the summaries to different txt files automatically. instead of copy-pasting, you can just run the following command for each model mod. the output will be saved to the file capturouputput.txt

# ----------------------------------------------------------------------------------------

#------------------------

###OUTPUT
# A DESeqDataSet is returned which contains all the fitted information within it, and the following section describes how to extract out results tables of interest from this object. This function will print out a message for the various steps it performs.  These are described in more detail in the manual page for DESeq, which can be accessed by typing ?DESeq. Briefly these are: the estimation of size factors (which control for differences in the library size of the sequencing drkeriments), the estimation of dispersion for each gene, and fitting a generalized linear model.



# DESeq2 normalised data (DESeq2 normalisation). 
# If we divide each column of the count table by the size factor for this column, 
# the count values are brought to a common scale, making them t8arable. When called with 
# normalized=TRUE, the counts accessor function [counts( cds, normalized=TRUE )] 
# performs this calculation. This is useful, e.g., for visualization."
write.csv(counts(dds,normalized=TRUE),file="TS1015_int_dds_normalized_counts_DESeq2.csv")
#original count data you imported can be given back here by specifying normalised=FALSE
write.csv(counts(dds,normalized=FALSE),file="TS1015_int_dds_NotNormalized_counts_DESeq2.csv")




####################################################################################
# bioinformatic analysis of RNA-seq data using DESeq2
# Code by Tahsha E. Say
# Oct, 2016
#
# Part b) 
# Transformations of raw count data 
#
#
####################################################################################


# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
#
#
# rlog stabilization and variance stabiliazation
# 
# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------


# transform raw counts into normalized values
# DESeq2 has two options:  
# 1) rlog transformed and 
# 2) variance stabilization -  good for heatmaps, etc.


# The two functions, rlog and varianceStabilizingTransformation, have an argument blind ie. the transformation should be blind to the sample information specied by the design formula.     

# blind so that the initial times setting does not influence the outcome, ie we want to see if the times cluster based purely on the individual datasets, in an unbiased way.  blind=TRUE should be used for t8aring samples in an manner unbiased by prior information on samples, for example to perform sample QA (quality assurance). blind=FALSE should be used for transforming data for downstream analysis, where the full use of the design information should be made. If many of genes have large differences in counts due to the drkerimental design, it is important to set blind=FALSE for downstream analysis.


# clear workspace - manually
#getwd() # get working directory
# setwd("~/path/to/working/directory/")
# directory <- "/path/to/counts/directory/"
setwd("/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_transformed_data") 



#########################################################################################
#-----------------------

# Can run prev scripts or 

# Alternatively read in output from step a03.03 (saved on 2016.10.26)

load(file = "/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_R_Scripts_wd_03.00/TS1015_h03.04_DESeq2_DEA_Design-group_2016.10.31/Environment_TS1015_h03.03_DEA_setup_Design-group_2016.10.31.RData")
#may need to load a aswell for design 
#saved in wd here

#-----------------------

# Now we want to transform the raw discretely distributed counts so that we can do clustering. To avoid that the distance measure is dominated by a few highly variable genes, and have a roughly equal contribution from all genes, we use it on the rlog-transformed data

rld <- rlogTransformation(dds, blind=TRUE) 
rlogMat <- assay(rld)
head(rlogMat)
#write.csv(rlogMat,file="TS1015_DESeq2_rlogMat_Blind_2016.11.21.csv")
write.table(rlogMat, file="TS1015_DESeq2_rlogMat_Blind_2016.11.21.txt", quote=FALSE, sep="\t")

# NOT blind
rldNTBL <- rlogTransformation(dds, blind=FALSE) 
rlogMatNTBL <- assay(rldNTBL)
head(rlogMatNTBL)
write.csv(rlogMatNTBL,file="TS1015_DESeq2_rlogMat_NotBlind_group_2016.11.21.csv")

#-----------------------

# transform the data without using any info specified by design formula (blind)
vsd<- varianceStabilizingTransformation(dds, blind=TRUE)
vsdMat <- assay(vsd)
#write.csv(vsdMat,file="TS1015_DESeq2_vsdMat_Blind_2016.11.21.csv")
write.table(vsdMat, file="TS1015_DESeq2_vsdMat_Blind_2016.11.21.txt", quote=FALSE, sep="\t")


# # # # transform the data USING info specified by design formula (NOT BLIND)
# NOT blind - use for heatmaps + comparision
vsdNTblind<- varianceStabilizingTransformation(dds, blind=FALSE)
vsdMatNTblind <- assay(vsdNTblind)
write.csv(vsdMatNTblind,file="TS1015_DESeq2_vsdMat_NotBlind_group_2016.11.21.csv")


save.image("/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_transformed_data/Environment_TS1015_b04.00_transfromations_group_2016.11.21.Rdata") ## to change to my wd




##########################################################################################
# 
# Reading in previously transformed data because the transformation step is time consuming 
# Name must include Mat (matrix) ie excel file. 
#
##########################################################################################

# comment these steps if NOT redoing analyses from this point (and comment the transformation step above). 
# must be in txt format (couldn't get the read.delim to work on 2016.04.13)
# rld = load(file = "./CEL-Seq-Sponge-networkConstruction-stepByStep.RData");





####################################################################################
# bioinformatic analysis of Cel-seq2 data using DESeq2
# Code by Tahsha E. Say
# Oct, 2016
#
# Part c)
#
#
####################################################################################



getwd() # get working directory

# setwd("~/path/to/working/directory/")
# directory <- "/path/to/counts/directory/"
setwd("/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_R_Scripts_wd_03.00/TS1015_c03.03_DESeq2_DEA_Design-LRT_2016.10.31") ## to change
getwd()


# notes
# Can run below script before starting: 
#/Users/tahshasay/Documents/Scripts/R_Scripts_TS1015/R_script_DESeq2_TS1015_03.00_2016/R_Script_TS1015_a03.04_DESeq2_DEA_setup_tm.trt.int_2016.10.28.sh


# Alternatively read in output from step a03.03 (saved on 2016.10.31)
load(file =  "/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_R_Scripts_wd_03.00/TS1015_a03.04_DEseq2_DEA_setup_Design-tm.trt.int_res_2016.10.31/Environment_TS1015_a03.04_DESeq2_DEA_setup_Design-tm.trt.int_2016.10.31.RData")


library("DESeq2") # start R session and load the DESeq2 package
library("ggplot2") # figures


# --------------------------------------------------------------------------------------

##########################################################################################
###
### gut of the DESeq2 analysis
###
##########################################################################################

# dds <- DESeq(dds)

# # # with no arguments to results, the results will be for the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable over the first level 
# ie. treatment: natural v. light
# he text time treated v untreated tells you the estimates are of log fold changed logs (trreated/ untreated). 
# res <- results(dds)


# This function will print out a message for the various steps it performs. 
# These are described in more detail in the manual page for DESeq, which can be accessed 
# by typing ?DESeq. Briefly these are: the estimation of size factors (which control for 
# differences in the library size of the sequencing drkeriments), 
# the estimation of dispersion for each gene, and fitting a generalized linear model
# A DESeqDataSet is returned which contains all the fitted information within it, 
# and the following section describes how to extract out results tables of interest from 
# this object.

# calling results without any argument will extract the estimated log2 fold changes and p 
# values for the last variable in the design formula. If there are more than two levels for 
# this variable, results will extract the results table for a t8arison of the last level 
# over the first level


# run the LRT analysis
dds <- DESeq(dds, test="LRT", reduced= ~1) ## to change


resultsNames(dds)
# output example: "Intercept" "Time_T45_vs_Tbase" "Time_Tend_vs_Tbase" "Genotype_mutant_vs_WT" "TimeT45.Genotypemutant" [6] "TimeTend.Genotypemutant"


TS1015_LRT_res <- results(dds) 

# with no arguments to results, the results will be for the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable over the first level. ie. treatment: natural v. light. The text time treated v untreated tells you the estimates are of log fold changed logs (trreated/ untreated). 
write.csv(as.data.frame(TS1015_LRT_res),file="TS1015_LRT_res_lncRNAs.csv")
dim(TS1015_LRT_res)
head(TS1015_LRT_res)
dim(TS1015_LRT_res)

# order results by padj value (most significant to least)
TS1015_Ordered <- TS1015_LRT_res[order(TS1015_LRT_res$padj),]
write.csv(as.data.frame(TS1015_Ordered),file="TS1015_LRT_2016.10.18_res_lncRNAs_resOrdered.csv")

dim(TS1015_Ordered)
head(TS1015_Ordered)
dim(TS1015_Ordered)



###OUTPUT
# A DESeqDataSet is returned which contains all the fitted information within it, and the following section describes how to extract out results tables of interest from this object. This function will print out a message for the various steps it performs.  These are described in more detail in the manual page for DESeq, which can be accessed by typing ?DESeq. Briefly these are: the estimation of size factors (which control for differences in the library size of the sequencing drkeriments), the estimation of dispersion for each gene, and fitting a generalized linear model.



# DESeq2 normalised data (DESeq2 normalisation). 
# If we divide each column of the count table by the size factor for this column, 
# the count values are brought to a common scale, making them t8arable. When called with 
# normalized=TRUE, the counts accessor function [counts( cds, normalized=TRUE )] 
# performs this calculation. This is useful, e.g., for visualization."
write.csv(counts(dds,normalized=TRUE),file="TS1015_LRT_dds_normalized_counts_DESeq2.csv")
#original count data you imported can be given back here by specifying normalised=FALSE
write.csv(counts(dds,normalized=FALSE),file="TS1015_LRT_dds_NotNormalized_counts_DESeq2.csv")



######################################################
## filter for sig genes with a 10% false discovery rate 
######################################################
# sig <- res[res$padj < 0.1,]
TS1015_Ordered_Sig_adjp0.1 <- subset(TS1015_Ordered, padj < 0.1)
write.csv(as.data.frame(TS1015_Ordered_Sig_adjp0.1),file="TS1015_LRT_2016.10.18_res_lncRNAs_resOrdered_Sig_adjp_0.1.csv")

write.table(as.data.frame(TS1015_Ordered_Sig_adjp0.1),file="TS1015_LRT_2016.10.18_res_lncRNAs_resOrdered_Sig_adjp_0.1.txt", quote=FALSE, sep="\t")

# check dimensions (no. diff LRT genes here) match in the excel sheet (no. should have)
dim(TS1015_Ordered_Sig_adjp0.1 )
mcols(TS1015_Ordered_Sig_adjp0.1 , use.names=TRUE)

######################################################
## filter for sig genes with a 5% false discovery rate
###################################################### 
# sig <- res[res$padj < 0.05,]
TS1015_Ordered_Sig_adjp0.05 <- subset(TS1015_Ordered, padj < 0.05)

write.csv(as.data.frame(TS1015_Ordered_Sig_adjp0.05),file="TS1015_LRT_2016.10.18_res_lncRNAs_resOrdered_Sig_adjp_0.05.csv")

write.table(as.data.frame(TS1015_Ordered_Sig_adjp0.05),file="TS1015_LRT_2016.10.18_res_lncRNAs_resOrdered_Sig_adjp_0.05.txt", quote=FALSE, sep="\t")

# check dimensions (no. diff LRT genes here) match in the excel sheet (no. should have)
dim(TS1015_Ordered_Sig_adjp0.05)
mcols(TS1015_Ordered_Sig_adjp0.05, use.names=TRUE)



# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------

########################################################
# print information needed for my excel sheet      # # # 
########################################################

resultsNames(dds)
dim(TS1015_LRT_res)
writeLines(capture.output(dim(TS1015_LRT_res)), "metadata_TS1015_LRT_res.txt")

dim(TS1015)
dim(expTS1015)

resultsNames(dds)

Design
sessionInfo()

# Capture the screen output into a character vector and use writeLines.
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

writeLines(capture.output(dim(TS1015_Ordered_Sig_adjp0.1)), "metadata_TS1015_LRT_Ordered_Sig_adjp0.1.txt")

writeLines(capture.output(dim(TS1015_Ordered_Sig_adjp0.05)), "metadata_TS1015_LRT_Ordered_Sig_adjp0.05.txt")


#‘sink’ diverts R output to a connection.
sink("sessionInfo_sink.txt")
sessionInfo()
sink()


# ---------------------------------------------------
# Save results for use in subsequent parts/ scripts
save(Design, dds, TS1015_LRT_res, file = "Environment_TS1015_c03.03_DEA_setup_LRT_2016.10.31.RData")
# --------------------------------------------------- date will change automatically

# to load above data in next script:
#__MEs= load(file = "./Environment_TS1015_a03.03_DEA_setup_tm.trt.int_2016.10.31.RData")

# ----------------------------------------------------------------------------------------





####################################################################################
# bioinformatic analysis of Cel-seq2 data using DESeq2
# Code by Tahsha E. Say
# Oct, 2016
#
# Part d) PCAs
#
#
####################################################################################



# add library dependencies
library("DESeq2") # start R session and load the DESeq2_2016.11.28 package
library("ggplot2")


setwd("/opsin/u/tahsha/TSay_R_wd/dXX.XX");


## to change
# read in previously transformed data (c04.06*) or uncomment steps below

# load above data in next script:
# read in previously transformed data (c04.06*) or uncomment steps below

# Alternatively read in output from step a03.03 (saved on 2016.10.26)
load(file = "Environment_TS1015_b03.01_transformations_2016.10.31.RData")
#load("/opsin/u/tahsha/TSay_R_wd/TS0616_ind_wNonInd_abc33.01_LRT/TS0616_*b33.01_LRT_ind_wNonInd_woBatch_vsd_2018.04.19.RData")
#vsd, rld vsdNTblind etc see a


#load(file = "/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_R_Scripts_wd_03.00/TS1015_q02.00_splsda_2017.06.14_tuning/splsda_2017.07.26_q03.01iii/Environment_TS1015_splsda_c_2017.06.19.RData") # to change


#vsd=read.table("/opsin/u/tahsha/TSay_R_wd/TS1015_transformed_data/TS1015_DESeq2_vsdMat_Blind_2016.11.21.txt", header=TRUE,row.names=1);



head(vsd)
#names(vsd)
colnames(vsd)
colData(vsd)


#----------------------------------------------------------
#basic PCA plot

#pdf("TS0616_ind_wNonInd_plotPCA_2018.03.13_time_NEW.pdf",8,6)
#plotPCA(vsd, "time")
#dev.off()
#----------------------------------------------------------


levels(vsd$time)
vsd$time <- factor(vsd$time, levels = c("3.5", "5", "7", "10"))


#col.list = c("sienna1", "firebrick1", "maroon4", "navy")
col.list = c("darkorange", "firebrick1", "magenta4", "blue")

# check variables
colData(vsd)


levels(vsd$treatment)


annot_col = Design

# col <- c("black","orange", "springgreen4") ## to change

ann_colours = list(
    treatment = c(lgt = "orange", nat ="springgreen4", drk = "black"),
    time = c(t1 = "lightgrey", t9 = "pink")
)


#-------------------------------------------------------------------------------------
pdf("FigX_TS1015_c03.03_PCA_plot_time_hpe_2018.04.27_grid_abbrev.pdf")
#print(plotPCA(vsd, intgroup=c("treatment")))	# this plots another pca in the same fileß	
# intgroup are the interesting groups for labeling the samples
# they tell the function to use them to choose colours

### remember to end the line with an operator so R knows that somepthing is coming (you are spreading the code over multiple lines)
  
data <- plotPCA(vsd, intgroup="time", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

ggplot(data, aes(PC1, PC2, color=group)) + 
geom_point(size=2.5) +
geom_point(aes(shape = vsd$treatment)) + #SHAPES
#geom_text(aes(label=colnames(vsd)), size=2.5, hjust=0.25, vjust=-0.6, show.legend = F) + # textual annotations http://docs.ggplot2.org/current/geom_text.html
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
stat_ellipse(type = "t") +
#ggtitle("Principal Components") +
#scale_color_brewer(palette="Dark2") +
#scale_colour_manual(values = col.list) +
scale_colour_manual(values = c("sienna1", "blue"),
			     name="Time (hours \npost emergence)",
                       breaks=c("t1", "t9"),
                       labels=c("1", "9")) +
scale_shape_discrete(name  ="Light regime",
                          breaks=c("nat", "drk", "lgt"),
                          labels=c("nat", "drk", "lgt")) +
#theme_classic() # White background, no gridlines
theme_bw() + # White background with grid lines
# theme_minimal() # Minimal theme - no border
# theme_grey() # grey background (default theme) 
#opts(panel.background = theme_rect(fill='white', colour='black')) + 
#opts(panel.grid.major = none, panel.grid.minor)
theme(axis.line = element_line(colour = "black"),
  #panel.grid.major.x = element_blank(),
  #panel.grid.major.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.background = element_blank(),
  legend.justification = "top",
)

dev.off() 



#-------------------------------------------------------------------------------------
pdf("FigXX_TS1015_c03.03_PCA_plot_time_hpe_2018.04.27_grid_abbrev.pdf")
#print(plotPCA(vsd, intgroup=c("treatment")))	# this plots another pca in the same fileß	
# intgroup are the interesting groups for labeling the samples
# they tell the function to use them to choose colours

### remember to end the line with an operator so R knows that somepthing is coming (you are spreading the code over multiple lines)
  
data <- plotPCA(vsd, intgroup="treatment", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

ggplot(data, aes(PC1, PC2, color=group)) + 
#geom_point(size=2.5) +
geom_point(aes(shape = vsd$time)) + #SHAPES
#geom_text(aes(label=colnames(vsd)), size=2.5, hjust=0.25, vjust=-0.6, show.legend = F) + # textual annotations http://docs.ggplot2.org/current/geom_text.html
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
stat_ellipse(type = "t") +
#ggtitle("Principal Components") +
#scale_color_brewer(palette="Dark2") +
#scale_colour_manual(values = col.list) +
scale_colour_manual(values = c("springgreen4", "orange", "black"),
							name ="Light regime",
                          breaks=c("nat", "drk", "lgt"),
                          labels=c("Natural day-night cycle", "Constant dark", "Constant light")) +
scale_shape_discrete(name="Time (hours \npost emergence)",
                       breaks=c("t1", "t9"),
                       labels=c("1", "9")) +
#theme_classic() # White background, no gridlines
theme_bw() + # White background with grid lines
# theme_minimal() # Minimal theme - no border
# theme_grey() # grey background (default theme) 
#opts(panel.background = theme_rect(fill='white', colour='black')) + 
#opts(panel.grid.major = none, panel.grid.minor)
theme(axis.line = element_line(colour = "black"),
  #panel.grid.major.x = element_blank(),
  #panel.grid.major.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.background = element_blank(),
  legend.justification = "top",
)

dev.off() 

# editing legends
# http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/
  #legend.title=element_blank() #hide the legend title / http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/

# legend.title=element_blank()




# natural day-night cycle
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
# select only nat treatment

selectNATURAL <-  c("r0t1nat1_01","r0t1nat2_02","r0t1nat3_03","r0t1nat4_04","r0t1nat5_05","r0t1nat6_06","r0t1nat7_07","r0t9nat1_22","r0t9nat2_23","r0t9nat3_24","r0t9nat4_25","r0t9nat5_43","r0t9nat6_26","r0t9nat7_27","r0t9nat8_28")
vsdnat <- vsd[,selectNATURAL]
rldnat <- rld[,selectNATURAL]
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------
pdf("FigX_nat_TS1015_c03.03_PCA_plot_time_hpe_2018.04.27_grid_abbrev.pdf")
#print(plotPCA(vsd, intgroup=c("treatment")))	# this plots another pca in the same fileß	
# intgroup are the interesting groups for labeling the samples
# they tell the function to use them to choose colours

### remember to end the line with an operator so R knows that somepthing is coming (you are spreading the code over multiple lines)
  
data <- plotPCA(vsdnat, intgroup="time", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

ggplot(data, aes(PC1, PC2, color=group)) + 
geom_point(size=2.5) +
#geom_point(aes(shape = vsdnat$treatment)) + #SHAPES
#geom_text(aes(label=colnames(vsd)), size=2.5, hjust=0.25, vjust=-0.6, show.legend = F) + # textual annotations http://docs.ggplot2.org/current/geom_text.html
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
stat_ellipse(type = "t") +
#ggtitle("Principal Components") +
#scale_color_brewer(palette="Dark2") +
#scale_colour_manual(values = col.list) +
scale_colour_manual(values = c("sienna1", "blue"),
			     name="Time (hours \npost emergence)",
                       breaks=c("t1", "t9"),
                       labels=c("1", "9")) +
#theme_classic() # White background, no gridlines
theme_bw() + # White background with grid lines
# theme_minimal() # Minimal theme - no border
# theme_grey() # grey background (default theme) 
#opts(panel.background = theme_rect(fill='white', colour='black')) + 
#opts(panel.grid.major = none, panel.grid.minor)
theme(axis.line = element_line(colour = "black"),
  #panel.grid.major.x = element_blank(),
  #panel.grid.major.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.background = element_blank(),
  legend.justification = "top",
)

dev.off() 

# editing legends
# http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/
  #legend.title=element_blank() #hide the legend title / http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/

# legend.title=element_blank()



# constant dark
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
# select only dark treatment
selectDARK <-  c("r0t1drk1_08","r0t1drk2_09","r0t1drk3_10","r0t1drk4_11","r0t1drk5_12","r0t1drk6_13","r0t1drk7_14","r0t9drk1_29","r0t9drk2_30","r0t9drk3_31","r0t9drk4_32","r0t9drk5_33","r0t9drk6_34","r0t9drk7_35")

rlddark <- rld[,selectDARK]
vsddark <- vsd[,selectDARK]
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
pdf("FigX_drk_TS1015_c03.03_PCA_plot_time_hpe_2018.04.27_grid_abbrev.pdf")
#print(plotPCA(vsd, intgroup=c("treatment")))	# this plots another pca in the same fileß	
# intgroup are the interesting groups for labeling the samples
# they tell the function to use them to choose colours

### remember to end the line with an operator so R knows that somepthing is coming (you are spreading the code over multiple lines)
  
data <- plotPCA(vsddrk, intgroup="time", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

ggplot(data, aes(PC1, PC2, color=group)) + 
geom_point(size=2.5) +
#geom_point(aes(shape = vsdnat$treatment)) + #SHAPES
#geom_text(aes(label=colnames(vsd)), size=2.5, hjust=0.25, vjust=-0.6, show.legend = F) + # textual annotations http://docs.ggplot2.org/current/geom_text.html
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
stat_ellipse(type = "t") +
#ggtitle("Principal Components") +
#scale_color_brewer(palette="Dark2") +
#scale_colour_manual(values = col.list) +
scale_colour_manual(values = c("sienna1", "blue"),
			     name="Time (hours \npost emergence)",
                       breaks=c("t1", "t9"),
                       labels=c("1", "9")) +
#theme_classic() # White background, no gridlines
theme_bw() + # White background with grid lines
# theme_minimal() # Minimal theme - no border
# theme_grey() # grey background (default theme) 
#opts(panel.background = theme_rect(fill='white', colour='black')) + 
#opts(panel.grid.major = none, panel.grid.minor)
theme(axis.line = element_line(colour = "black"),
  #panel.grid.major.x = element_blank(),
  #panel.grid.major.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.background = element_blank(),
  legend.justification = "top",
)

dev.off() 



# constant light
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
# select only lgt treatment

selectLIGHT <-  c("r0t1lgt1_15","r0t1lgt2_16","r0t1lgt3_17","r0t1lgt4_18","r0t1lgt5_19","r0t1lgt6_20","r0t1lgt7_21","r0t9lgt1_36","r0t9lgt2_37","r0t9lgt3_38","r0t9lgt4_39","r0t9lgt5_40","r0t9lgt6_41","r0t9lgt7_42")
vsdlgt <- vsd[,selectLIGHT]
rldlgt <- rld[,selectLIGHT]
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
pdf("FigX_lgt_TS1015_c03.03_PCA_plot_time_hpe_2018.04.27_grid_abbrev.pdf")
#print(plotPCA(vsd, intgroup=c("treatment")))	# this plots another pca in the same fileß	
# intgroup are the interesting groups for labeling the samples
# they tell the function to use them to choose colours

### remember to end the line with an operator so R knows that somepthing is coming (you are spreading the code over multiple lines)
  
data <- plotPCA(vsdlgt, intgroup="time", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

ggplot(data, aes(PC1, PC2, color=group)) + 
#geom_point(size=2.5) +
#geom_point(aes(shape = vsdnat$treatment)) + #SHAPES
#geom_text(aes(label=colnames(vsd)), size=2.5, hjust=0.25, vjust=-0.6, show.legend = F) + # textual annotations http://docs.ggplot2.org/current/geom_text.html
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
stat_ellipse(type = "t") +
#ggtitle("Principal Components") +
#scale_color_brewer(palette="Dark2") +
#scale_colour_manual(values = col.list) +
scale_colour_manual(values = c("sienna1", "blue"),
			     name="Time (hours \npost emergence)",
                       breaks=c("t1", "t9"),
                       labels=c("1", "9")) +
#theme_classic() # White background, no gridlines
theme_bw() + # White background with grid lines
# theme_minimal() # Minimal theme - no border
# theme_grey() # grey background (default theme) 
#opts(panel.background = theme_rect(fill='white', colour='black')) + 
#opts(panel.grid.major = none, panel.grid.minor)
theme(axis.line = element_line(colour = "black"),
  #panel.grid.major.x = element_blank(),
  #panel.grid.major.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.background = element_blank(),
  legend.justification = "top",
)

dev.off() 

# editing legends
# http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/
  #legend.title=element_blank() #hide the legend title / http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/

# legend.title=element_blank()
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


########################################################
# customise the PCA plot using the ggplot function # # # 
########################################################

pdf("TS0616_ind_wNonInd_d33.01_LRT_PCA_plot_vsd_trialshapes.pdf")
#print(plotPCA(vsd, intgroup=c("treatment"))) # this plots another pca in same file		
# intgroup are the interesting groups for labeling the samples
# they tell the function to use them to choose colours

### remember to end the line with an operator so R knows that somepthing is coming (you are spreading the code over multiple lines)
   
#-------------------------------------------------------
#http://ggplot2.tidyverse.org/reference/geom_point.html  
#p + geom_point(aes(shape = factor(cyl)))
#------------------------------------------------------- 

data <- plotPCA(vsd, intgroup=c("time"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

ggplot(data, aes(PC1, PC2, color=group)) + 
#geom_point(size=3) +
geom_point(aes(shape = vsd$treatment)) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
scale_colour_manual(values = col.list) +
#stat_ellipse(type = "t") +
ggtitle("Principal Components") +
#scale_color_brewer(palette="Dark2") +
#scale_colour_manual(values = col.list) +
#theme_classic() # White background, no gridlines
theme_bw() + # White background with grid lines
# theme_minimal() # Minimal theme - no border
# theme_grey() # grey background (default theme) 
#opts(panel.background = theme_rect(fill='white', colour='black')) + 
#opts(panel.grid.major = none, panel.grid.minor)
theme(axis.line = element_line(colour = "black"),
  panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.background = element_blank()
)

dev.off() 
          
          
########################################################
# w elipses  
########################################################

# Ellipses (ggplot2)
# http://ggplot2.tidyverse.org/reference/stat_ellipse.html
# The default "t" assumes a multivariate t-distribution


pdf("TS0616_ind_wNonInd_d33.01_LRT_PCA_plot_vsd_elipse.pdf") 
data <- plotPCA(vsd, intgroup=c("time"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

ggplot(data, aes(PC1, PC2, color=group)) + 
#geom_point(size=3) +
geom_point(aes(shape = vsd$treatment)) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
scale_colour_manual(values = col.list) +
stat_ellipse(type = "t") +
ggtitle("Principal Components") +
#scale_color_brewer(palette="Dark2") +
#scale_colour_manual(values = col.list) +
#theme_classic() # White background, no gridlines
theme_bw() + # White background with grid lines
# theme_minimal() # Minimal theme - no border
# theme_grey() # grey background (default theme) 
#opts(panel.background = theme_rect(fill='white', colour='black')) + 
#opts(panel.grid.major = none, panel.grid.minor)
theme(axis.line = element_line(colour = "black"),
  panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.background = element_blank()
)

dev.off() 



########################################################
# w labelled points  
########################################################

# textual annotations
#geom_point() +
#geom_text(aes(label=names), size=2.5, hjust=0.25, vjust=-0.6, show.legend = F) + # textual annotations http://docs.ggplot2.org/current/geom_text.html


pdf("TS0616_ind_wNonInd_d33.01_LRT_PCA_plot_vsd_elipse_labels.pdf")

  
data <- plotPCA(vsd, intgroup=c("time"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

ggplot(data, aes(PC1, PC2, color=group)) + 
#geom_point(size=3) +
geom_text(aes(label=colnames(vsd)), size=2.5, hjust=0.25, vjust=-0.6, show.legend = F) + # textual annotations http://docs.ggplot2.org/current/geom_text.html
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
scale_colour_manual(values = col.list) +
stat_ellipse(type = "t") +
ggtitle("Principal Components") +
#scale_color_brewer(palette="Dark2") +
#scale_colour_manual(values = col.list) +
#theme_classic() # White background, no gridlines
theme_bw() + # White background with grid lines
# theme_minimal() # Minimal theme - no border
# theme_grey() # grey background (default theme) 
#opts(panel.background = theme_rect(fill='white', colour='black')) + 
#opts(panel.grid.major = none, panel.grid.minor)
theme(axis.line = element_line(colour = "black"),
  panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.background = element_blank()
)

dev.off() 





####################################################################################
# bioinformatic analysis of Cel-seq2 data using DESeq2
# Code by Tahsha E. Say
# Oct, 2016
#
# 
# Part e) Heatmaps of DEG (previously named f)
#
####################################################################################


#------------------------------------------------------------
# Run below every time you want to run this script on piwi
#------------------------------------------------------------
# load wd
#
# ssh server
#
# cd ./TSay_R_wd/f07.00_LRT_heatmaps/
#  module avail|less |grep R
# 
# module load R/3.2.3-foss-2016a
#
#
# to run this self contained script
# nohup Rscript #R_Script_TS1015_f07.00c_LRT_DESeq2_heatmaps_DEG_2016.10.31_colour_2017.06.16_piwi.R &
#------------------------------------------------------------


setwd("/opsin/u/tahsha/TSay_R_wd/f07.00_LRT_heatmaps")
getwd()


# move output to:
#setwd("/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_R_Scripts_wd_03.00/TS1015_f03.05c_DESeq2_heatmaps_Design-LRT_2017.06.16");

getwd()


# add library dependencies
library("DESeq2") # start R session and load the DESeq2_2016.11.28 package
library("ggplot2")

# Alternatively read in output from step a03.03 (saved on 2016.10.26)
load(file = "Environment_TS1015_c03.03_DEA_setup_LRT_2016.10.31.RData")



# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------

# HEATMAPS OF DIFFERENTIALLY EXPRESSED GENES on transformed data (BLIND)

# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------


# comment these steps if NOT redoing analyses from this point 

###################
# read in vsd data - to change
###################
# comment these if not starting analyses here. 
#vsd=read.table("/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_transformed_data/TS1015_DESeq2_group_vsdMat_Blind_2016.10.24_Rv3.3.1_inclHeader.txt", header=TRUE,row.names=1);

vsd=read.table("../TS1015_transformed_data/TS1015_DESeq2_vsdMat_Blind_2016.11.21.txt", header=TRUE,row.names=1);

###################
# read in rld data - to change
###################
#rld=read.table("/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_transformed_data/TS1015_DESeq2_group_rlogMat_Blind_2016.10.24_Rv3.3.1_inclHeader.txt", header=TRUE,row.names=1);

rld=read.table("../TS1015_transformed_data/TS1015_DESeq2_rlogMat_Blind_2016.11.21.txt", header=TRUE,row.names=1);

#vsdNTBL=read.table("/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_transformed_data/TS1015_DESeq2_group_vsdMat_NotBlind_2016.10.24_Rv3.3.1.txt", header=TRUE,row.names=1)

#vsdNTBL=read.table("TS1015_DESeq2_vsdMat_NotBlind_group_2016.11.21.txt", header=TRUE,row.names=1)

# Alternatively
# rld = load(file = "./CEL-Seq-Sponge-networkConstruction-stepByStep.RData");

#----------------------------------------------------


##############################
# re-order TRANSFORMED DATA used to create heatmaps to make ordering intuitive.
# best to have "nat" treatment in the middle (close to either light or dark dep on time)
# Here I use R to Reorder matrix columns by matching colnames to list of string
# Ref: http://stackoverflow.com/questions/25446714/r-reorder-matrix-columns-by-matching-colnames-to-list-of-string

col.order <- c("r0t1lgt1_15","r0t1lgt2_16","r0t1lgt3_17","r0t1lgt4_18","r0t1lgt5_19","r0t1lgt6_20","r0t1lgt7_21","r0t1nat1_01","r0t1nat2_02","r0t1nat3_03","r0t1nat4_04","r0t1nat5_05","r0t1nat6_06","r0t1nat7_07","r0t1drk1_08","r0t1drk2_09","r0t1drk3_10","r0t1drk4_11","r0t1drk5_12","r0t1drk6_13","r0t1drk7_14","r0t9lgt1_36","r0t9lgt2_37","r0t9lgt3_38","r0t9lgt4_39","r0t9lgt5_40","r0t9lgt6_41","r0t9lgt7_42","r0t9nat1_22","r0t9nat2_23","r0t9nat3_24","r0t9nat4_25","r0t9nat5_43","r0t9nat6_26","r0t9nat7_27","r0t9nat8_28","r0t9drk1_29","r0t9drk2_30","r0t9drk3_31","r0t9drk4_32","r0t9drk5_33","r0t9drk6_34","r0t9drk7_35")

vsd <- vsd[,col.order]
rld <- vsd[,col.order]



##########################################################################################
#
# Choose heatmap colours
#
##########################################################################################

library(ggplot2)
library(gplots) 
library(pheatmap)
library(RColorBrewer)
library("colorRamps") # for blue2red 


### remember to scale heat map by row (especially if you sorted by total expression level).  This will make the trend (down over different genes for each larva less obvious).  But will make your Left to Right trend MORE obvious.


##############################
# Diverging colours
##############################

#hmcol <- colorRampPalette(c("blue",'white','firebrick2'))(n=1000)
hmcol = blue2red(400) # deceptive becuase lightt blue and yellow are sim exp level.
#hmcol <- colorRampPalette(rev(brewer.pal(9, "Spectral")))(1000)
# default used for heatmpa3.  
#hmcol <- colorRampPalette(c("blue",'white','red'))(n=1000)


##############################
# Sequential
##############################
#hmcol <- colorRampPalette(brewer.pal(9, "Blues")) (2000)

#---------------
#library(viridis)
#hmcol <- rev(viridis(256))
#hmcol <- magma(256)
#hmcol <- rev(plasma(256))
#inferno
#---------------



annot_col = Design

# col <- c("black","orange", "springgreen4") ## to change

ann_colours = list(
    treatment = c(lgt = "orange", nat ="springgreen4", drk = "black"),
    time = c(t1 = "lightgrey", t9 = "pink")
)



##########################################################################################
#
# Heatmap DEG. Visualising differential expression using heat maps
# 
# all DEG (adjp<0.1)
#
##########################################################################################

# transform data without using any info specified by design formula (BLIND)


selectDEG <- order(TS1015_LRT_res$padj)[1:7061] ## to change (no. DEG)

########################################################
# vsd    
########################################################

pdf("TS1015_2016.10.26_LRT_DEG7061_adjp0.1_vsd_scalerow_clusteringboth.pdf", width = 8, height = 10,  onefile=FALSE)
#colnames(vsd) <- paste( rld$time, 1:23, sep="-" )
pheatmap(vsd[selectDEG,], col=hmcol, 
cluster_cols=TRUE, 
cluster_rows=TRUE, 
show_rownames=FALSE, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off()

#cluster cols
pdf("TS1015_2016.10.26_LRT_DEG7061_adjp0.1_vsd_scalerow_clusteringcols.pdf", width = 8, height = 10,  onefile=FALSE)
pheatmap(vsd[selectDEG,], col=hmcol, 
cluster_cols=TRUE, 
cluster_rows=FALSE, 
show_rownames=FALSE, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off() 

# cluster rows
pdf("TS1015_2016.10.26_LRT_DEG7061_adjp0.1_vsd_scalerow_clusteringrows.pdf", width = 8, height = 10,  onefile=FALSE)
pheatmap(vsd[selectDEG,], col=hmcol, 
cluster_cols=FALSE, 
cluster_rows=TRUE, 
show_rownames=FALSE, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off()

########################################################
# rld    
########################################################

pdf("TS1015_2016.10.26_LRT_DEG7061_adjp0.1_rld_scalerow_clusteringboth.pdf", width = 8, height = 10,  onefile=FALSE)
#colnames(vsd) <- paste( rld$time, 1:23, sep="-" )
pheatmap(rld[selectDEG,], col=hmcol, 
cluster_cols=TRUE, 
cluster_rows=TRUE, 
show_rownames=FALSE, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off()

#cluster cols
pdf("TS1015_2016.10.26_LRT_DEG7061_adjp0.1_rld_scalerow_clusteringcols.pdf", width = 8, height = 10,  onefile=FALSE)
pheatmap(rld[selectDEG,], col=hmcol, 
cluster_cols=TRUE, 
cluster_rows=FALSE, 
show_rownames=FALSE, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off() 

# cluster rows
pdf("TS1015_2016.10.26_LRT_DEG7061_adjp0.1_rld_scalerow_clusteringrows.pdf", width = 8, height = 10,  onefile=FALSE)
pheatmap(rld[selectDEG,], col=hmcol, 
cluster_cols=FALSE, 
cluster_rows=TRUE, 
show_rownames=FALSE, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off()



##########################################################################################
#
# Heatmap DEG. Visualising differential expression using heat maps
# 
# all DEG (adjp<0.05)
#
##########################################################################################

# transform data without using any info specified by design formula (BLIND)
#
#
#

selectDEG <- order(TS1015_LRT_res$padj)[1:5737] ## to change (no. DEG) this is all at the adjp<0.05 cut off


########################################################
# vsd    
########################################################

pdf("TS1015_2016.10.26_LRT_DEG5737_adjp0.05_vsd_scalerow_clusteringboth.pdf", width = 8, height = 10,  onefile=FALSE)
#colnames(vsd) <- paste( rld$time, 1:23, sep="-" )
pheatmap(vsd[selectDEG,], col=hmcol, 
cluster_cols=TRUE, 
cluster_rows=TRUE, 
show_rownames=FALSE, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off()

#cluster cols
pdf("TS1015_2016.10.26_LRT_DEG5737_adjp0.05_vsd_scalerow_clusteringcols.pdf", width = 8, height = 10,  onefile=FALSE)
pheatmap(vsd[selectDEG,], col=hmcol, 
cluster_cols=TRUE, 
cluster_rows=FALSE, 
show_rownames=FALSE, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off() 

# cluster rows
pdf("TS1015_2016.10.26_LRT_DEG5737_adjp0.05_vsd_scalerow_clusteringrows.pdf", width = 8, height = 10,  onefile=FALSE)
pheatmap(vsd[selectDEG,], col=hmcol, 
cluster_cols=FALSE, 
cluster_rows=TRUE, 
show_rownames=FALSE, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off()

########################################################
# rld    
########################################################

pdf("TS1015_2016.10.26_LRT_DEG5737_adjp0.05_rld_scalerow_clusteringboth.pdf", width = 8, height = 10,  onefile=FALSE)
#colnames(vsd) <- paste( rld$time, 1:23, sep="-" )
pheatmap(rld[selectDEG,], col=hmcol, 
cluster_cols=TRUE, 
cluster_rows=TRUE, 
show_rownames=FALSE, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off()

#cluster cols
pdf("TS1015_2016.10.26_LRT_DEG5737_adjp0.05_rld_scalerow_clusteringcols.pdf", width = 8, height = 10,  onefile=FALSE)
pheatmap(rld[selectDEG,], col=hmcol, 
cluster_cols=TRUE, 
cluster_rows=FALSE, 
show_rownames=FALSE, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off() 


# cluster rows
pdf("TS1015_2016.10.26_LRT_DEG5737_adjp0.05_rld_scalerow_clusteringrows.pdf", width = 8, height = 10,  onefile=FALSE)
pheatmap(rld[selectDEG,], col=hmcol, 
cluster_cols=FALSE, 
cluster_rows=TRUE, 
show_rownames=FALSE, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off()



#------------------------------------
# cut the tree: 2 
#------------------------------------

# cluster rows
pdf("TS1015_2016.10.26_LRT_DEG5737_adjp0.05_vsd_scalerow_clusteringrows_cut_pearsons.pdf", width = 8, height = 10,  onefile=FALSE)
pheatmap(vsd[selectDEG,], col=hmcol, 
cluster_cols=FALSE, 
cluster_rows=TRUE, 
clustering_distance_rows = "correlation", #Pearsons
show_rownames=FALSE, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
cutree_rows=2, #number of clusters the rows are divided into, based on the hierarchical clustering (using cutree), if rows are not clustered, the argument is ignored
scale="row")
dev.off()


# cluster rows
pdf("TS1015_2016.10.26_LRT_DEG5737_adjp0.05_vsd_scalerow_clusteringrows_cut_kmeans.pdf", width = 8, height = 10,  onefile=FALSE)
pheatmap(vsd[selectDEG,], col=hmcol, 
cluster_cols=FALSE, 
cluster_rows=TRUE, 
show_rownames=FALSE, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
kmeans_k = 2, # the number of kmeans clusters to make, if we want to agggregate the rows before
drawing heatmap. If NA then the rows are not aggregated
scale="row")
dev.off()

#------------------------------------
# cut the tree: 2 
#------------------------------------


res <- pheatmap(vsd, col=hmcol, 
cluster_cols=FALSE, 
cluster_rows=TRUE, 
show_rownames=FALSE, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")


cutree(res$tree_row, k = 2)

# returns the cluster membership for each row of your original data
clust <- cbind(vsd, cluster = cutree(res$tree_row, k = 2))

write.table(as.data.frame(clust),file="TS1015_2016.11.09_VENNY_DEG1248_vsd_scalerow_clusteringrows_shared.nat.drk_cut2.txt", quote=FALSE, sep="\t", col.names=NA)
# cluster = col of interest


#----------------------------------------

# To see it in the order as it is plotted in the heatmap, you have order it according to the heat map index. 


# grab a cluster
cluster1 <- clust[clust$cluster == 1,]

dim(cluster1)
# = 44 (and I want to remove the last "cluster" col cause I do not want this in my  heatmap).

write.table(as.data.frame(cluster1[,1:43]),file="TS1015_2016.11.09_VENNY_DEG1248_vsd_scalerow_clusteringrows_shared.nat.drk_cluster1.txt", quote=FALSE, sep="\t", col.names=NA)


pdf("TS1015_2016.11.09_VENNY_DEG1248_vsd_scalerow_clusteringrows_shared.nat.drk_cluster1.pdf", width = 8, height = 10,  onefile=FALSE)
pheatmap(cluster1[,1:43], col=hmcol, 
cluster_cols=FALSE, 
cluster_rows=TRUE, 
show_rownames=FALSE, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off()

#---------------------- write just Aq list
write.table(as.data.frame(rownames(cluster1)),file="TS1015_2016.11.09_VENNY_DEG1248_vsd_scalerow_clusteringrows_shared.nat.drk_cluster1_Aqlist.txt", quote=FALSE, eol = "\r", col.names=FALSE, row.names=FALSE)
#------------------------------------------



########################################################
# print information needed for my excel sheet      # # # 
########################################################


# Capture the screen output into a character vector and use writeLines.
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")


#‘sink’ diverts R output to a connection.
sink("sessionInfo_sink.txt")
sessionInfo()
sink()



####################################################################################
# bioinformatic analysis of RNA-seq data using DESeq2
# Code by Tahsha E. Say
# Oct, 2016
#
# Part f)
#
#
####################################################################################


setwd("/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_R_Scripts_wd_03.00/TS1015_h03.04_DESeq2_DEA_Design-group_2016.10.31") ## to change
getwd()

# Alternatively read in output from step a03.03 (saved on 2016.10.31)
load(file =  "/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_R_Scripts_wd_03.00/TS1015_a03.04_DEseq2_DEA_setup_Design-tm.trt.int_res_2016.10.31/Environment_TS1015_a03.04_DESeq2_DEA_setup_Design-tm.trt.int_2016.10.31.RData")

# Transformations were conducted on a diff day/ prev
# R_Script_TS1015_b03.00_DESeq2_transformations_lncRNA_2016.10.20.sh ** changed to todays date

# add library dependencies
library("DESeq2") # start R session and load the DESeq2 package
library("ggplot2") # figures



# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------
#
#  DEA: Pairwise comparisons - create new dds dataset!!! ~group !!!
#
# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------

 # Many users begin to add interaction terms to the design formula, when in fact a much simpler approach would give all the results tables that are desired. We will explain this approach first, because it is much simpler to perform. If the comparisons of interest are, for example, the effect of a condition for different sets of samples, a simpler approach than adding interaction terms explicitly to the design formula is to perform the following steps:
 
 # 1. combine the factors of interest into a single factor with all combinations of the original factors
 # 2. change the design to include just this factor, e.g.  group

# Using this design is similar to adding an interaction term, in that it models multiple condition effects which can be easily extracted with results. Suppose we have two factors genotype (with values I, II, and III) and condition (with values A and B), and we want to extract the condition effect specifically for each genotype. We could use the following approach to obtain, e.g. the condition effect for genotype I:


# create a new design formula 

dds$group <- factor(paste0(dds$time, dds$treatment))
design(dds) <- ~ group
dds <- DESeq(dds)
resultsNames(dds) # display the levels of the new factor you created


# DESeq2 normalised data (DESeq2 normalisation). 
# If we divide each column of the count table by the size factor for this column, 
# the count values are brought to a common scale, making them t8arable. When called with 
# normalized=TRUE, the counts accessor function [counts( cds, normalized=TRUE )] 
# performs this calculation. This is useful, e.g., for visualization."
write.csv(counts(dds,normalized=TRUE),file="TS1015_group_dds_normalized_counts_DESeq2.csv")

# Example save R data save.image("~/Documents/LABO/DEGNAN/PROYECTOS/GRN/DATA/WGCNA/CEL-SEQ-GRN/Environment_GRN_11052015.RData")   ##### CARMEL TO CHANGE
#original count data you imported can be given back here by specifying normalised=FALSE
write.csv(counts(dds,normalized=FALSE),file="TS1015_group_dds_NotNormalized_counts_DESeq2.csv")


# t1 v t9 (i.e. 1 hpe v 9 hpe)

t1.t9.drk <- results(dds, contrast=c("group", "t1drk", "t9drk")) 
t1.t9.lgt <- results(dds, contrast=c("group", "t1lgt", "t9lgt")) 
t1.t9.nat <- results(dds, contrast=c("group", "t1nat", "t9nat")) 

# summarise some basic tallies suing the summary function
summary(t1.t9.drk)
summary(t1.t9.lgt)
summary(t1.t9.nat)

# We can save the table, and also print out some information on what the columns mean
mcols(t1.t9.drk,use.names=TRUE)
mcols(t1.t9.lgt,use.names=TRUE)
mcols(t1.t9.nat,use.names=TRUE)

# more info about which variables and tests were used can be found by calling the function mcols on the res object
mcols(t1.t9.drk)$description
mcols(t1.t9.lgt)$description
mcols(t1.t9.nat)$description


write.csv(as.data.frame(mcols(t1.t9.drk,use.names=TRUE)),file = "metadata_cols_t1.t9.drk_2016.10.31.csv") ## to change

write.csv(as.data.frame(mcols(t1.t9.lgt, use.name = T)),file = "metadata_cols_t1.t9.lgt_2016.10.31.csv") ## to change

write.csv(as.data.frame(mcols(t1.t9.nat, use.name = T)),file = "metadata_cols_t1.t9.nat_2016.10.31.csv") ## to change


library(dplyr)
library(magrittr)

# order samples
#t1.t9.lgt_Ordered <- t1.t9.lgt[order(t1.t9.lgt$padj),]

# Get the results for the factorA vs factorB samples and arrange by p-value.
t1.t9.drk_Ordered <- results(dds, tidy=TRUE, contrast=c("group", "t1drk", "t9drk")) %>%
	arrange(padj, pvalue) %>%
	tbl_df() 
t1.t9.drk_Ordered

t1.t9.lgt_Ordered <- results(dds, tidy=TRUE, contrast=c("group", "t1lgt", "t9lgt")) %>%
	arrange(padj, pvalue) %>%
	tbl_df() 
t1.t9.lgt_Ordered

t1.t9.nat_Ordered <- results(dds, tidy=TRUE, contrast=c("group", "t1nat", "t9nat")) %>%
	arrange(padj, pvalue) %>%
	tbl_df() 
t1.t9.nat_Ordered


#-----------------------------------------------------------
# Deprecated - this is incorp above. 
# order results by padj value (most significant to least)
#t1.t9.drk_Ordered <- t1.t9.drk[order(t1.t9.drk$padj),]

#head(t1.t9.drk_Ordered)
write.csv(as.data.frame(t1.t9.drk_Ordered),file="TS1015_t1.t9.drk_resOrdered.csv")
write.table(as.data.frame(t1.t9.drk_Ordered),file="TS1015_t1.t9.drk_resOrdered.txt", quote=FALSE, sep="\t")
dim(t1.t9.drk_Ordered)
# 28061, 6

#t1.t9.lgt_Ordered <- t1.t9.lgt[order(t1.t9.lgt$padj),]
head(t1.t9.lgt_Ordered)
write.csv(as.data.frame(t1.t9.lgt_Ordered),file="TS1015_t1.t9.lgt_resOrdered.csv")
write.table(as.data.frame(t1.t9.lgt_Ordered),file="TS1015_t1.t9.lgt_resOrdered.txt", quote=FALSE, sep="\t")
dim(t1.t9.lgt_Ordered)
# 28061, 6

#t1.t9.nat_Ordered <- t1.t9.nat[order(t1.t9.nat$padj),]
head(t1.t9.nat_Ordered)
write.csv(as.data.frame(t1.t9.nat_Ordered),file="TS1015_t1.t9.nat_resOrdered.csv")
write.table(as.data.frame(t1.t9.nat_Ordered),file="TS1015_t1.t9.nat_resOrdered.txt", quote=FALSE, sep="\t")
dim(t1.t9.nat_Ordered)
# 28061, 6

#-----------------------------------------------------------
## filter for sig genes with a 10% false discovery rate 
# sig <- res[res$padj < 0.1,]
t1.t9.drk_OrderedSig <- subset(t1.t9.drk_Ordered, padj < 0.1)
write.csv(as.data.frame(t1.t9.drk_OrderedSig),file="TS1015_t1.t9.drk_resOrdered_Sig_adjp_0.1.csv")
write.table(as.data.frame(t1.t9.drk_OrderedSig),file="TS1015_t1.t9.drk_resOrdered_Sig_adjp_0.1.txt", quote=FALSE, sep="\t")

# check dimensions (no. diff drk genes here) match in the excel sheet (no. should have)
dim(t1.t9.drk_OrderedSig)


## filter for sig genes with a 5% false discovery rate 
# sig <- res[res$padj < 0.05,]
t1.t9.drk_OrderedSig_padj0.05 <- subset(t1.t9.drk_Ordered, padj < 0.05)
write.csv(as.data.frame(t1.t9.drk_OrderedSig_padj0.05),file="TS1015_t1.t9.drk_resOrdered_Sig_adjp_0.05.csv")
write.table(as.data.frame(t1.t9.drk_OrderedSig_padj0.05),file="TS1015_t1.t9.drk_resOrdered_Sig_adjp_0.05.txt", quote=FALSE, sep="\t")

# check dimensions (no. diff drk genes here) match in the excel sheet (no. should have)
dim(t1.t9.drk_OrderedSig_padj0.05)
mcols(t1.t9.drk_OrderedSig_padj0.05, use.names=TRUE)


#-----------------------------------------------------------
## filter for sig genes with a 10% false discovery rate 
# sig <- res[res$padj < 0.1,]
t1.t9.lgt_OrderedSig <- subset(t1.t9.lgt_Ordered, padj < 0.1)
write.csv(as.data.frame(t1.t9.lgt_OrderedSig),file="TS1015_t1.t9.lgt_resOrdered_Sig_adjp_0.1.csv")
write.table(as.data.frame(t1.t9.lgt_OrderedSig),file="TS1015_t1.t9.lgt_resOrdered_Sig_adjp_0.1.txt", quote=FALSE, sep="\t")


## filter for sig genes with a 5% false discovery rate 
# sig <- res[res$padj < 0.05,]
t1.t9.lgt_OrderedSig_padj0.05 <- subset(t1.t9.lgt_Ordered, padj < 0.05)
write.csv(as.data.frame(t1.t9.lgt_OrderedSig_padj0.05),file="TS1015_t1.t9.lgt_resOrdered_Sig_adjp_0.05.csv")
write.table(as.data.frame(t1.t9.lgt_OrderedSig_padj0.05),file="TS1015_t1.t9.lgt_resOrdered_Sig_adjp_0.05.txt", quote=FALSE, sep="\t")


#-----------------------------------------------------------

## filter for sig genes with a 10% false discovery rate 
# sig <- res[res$padj < 0.1,]
t1.t9.nat_OrderedSig <- subset(t1.t9.nat_Ordered, padj < 0.1)
write.csv(as.data.frame(t1.t9.nat_OrderedSig),file="TS1015_t1.t9.nat_resOrdered_Sig_adjp_0.1.csv")
write.table(as.data.frame(t1.t9.nat_OrderedSig),file="TS1015_t1.t9.nat_resOrdered_Sig_adjp_0.1.txt", quote=FALSE, sep="\t")



## filter for sig genes with a 5% false discovery rate 
# sig <- res[res$padj < 0.05,]
t1.t9.nat_OrderedSig_padj0.05 <- subset(t1.t9.nat_Ordered, padj < 0.05)
write.csv(as.data.frame(t1.t9.nat_OrderedSig_padj0.05),file="TS1015_t1.t9.nat_resOrdered_Sig_adjp_0.05.csv")
write.table(as.data.frame(t1.t9.nat_OrderedSig_padj0.05),file="TS1015_t1.t9.nat_resOrdered_Sig_adjp_0.05.txt", quote=FALSE, sep="\t")


########################################################
# print information needed for my excel sheet      # # # 
########################################################

resultsNames(dds)

Design
sessionInfo()

# Capture the screen output into a character vector and use writeLines.
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

writeLines(capture.output(dim(t1.t9.drk_OrderedSig)), "metadata_TS1015_t1.t9.drk_OrderedSig_adjp0.1.txt")

writeLines(capture.output(dim(t1.t9.nat_OrderedSig)), "metadata_TS1015_t1.t9.nat_OrderedSig_adjp0.1.txt")

writeLines(capture.output(dim(t1.t9.lgt_OrderedSig)), "metadata_TS1015_t1.t9.lgt_OrderedSig_adjp0.1.txt")

######

writeLines(capture.output(dim(t1.t9.drk_OrderedSig_padj0.05)), "metadata_TS1015_t1.t9.drk_OrderedSig_adjp0.05.txt")

writeLines(capture.output(dim(t1.t9.nat_OrderedSig_padj0.05)), "metadata_TS1015_t1.t9.nat_OrderedSig_adjp0.05.txt")

writeLines(capture.output(dim(t1.t9.lgt_OrderedSig_padj0.05)), "metadata_TS1015_t1.t9.lgt_OrderedSig_adjp0.05.txt")


#‘sink’ diverts R output to a connection.
sink("sessionInfo_sink.txt")
sessionInfo()
sink()


# ---------------------------------------------------
# Save results for use in subsequent parts/ scripts
save(dds, Design, t1.t9.drk, t1.t9.lgt, t1.t9.nat, t1.t9.drk_Ordered, t1.t9.nat_Ordered, t1.t9.lgt_Ordered,  file = "Environment_TS1015_h03.03_DEA_setup_group_2016.10.31.RData")
# --------------------------------------------------- date will change automatically

# load above file
#__MEs= load(file = "./Environment_TS1015_a03.03_DEA_setup_tm.trt.int_2016.10.31.RData")





####################################################################################
# bioinformatic analysis of Cel-seq2 data using DESeq2
# Code by Tahsha E. Say
# Oct, 2016
#
# 
# Part g) plot counts
#
#
# DESeq2 post analysis plots of expression values with a multipage PDF catalog	#
# Code by: Stephen Turner										#
# June 6, 2016												#												#
# https://rpubs.com/turnersd/plot-deseq-results-multipage-pdf				#
#
####################################################################################


# Goal: "make a boxplot of the expression values of a gene of interest"

getwd() # get working directory
#
#
# Create new wd for todays results. 
# setwd("~/path/to/working/directory/")
# directory <- "/path/to/counts/directory/"
setwd("/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_R_Scripts_wd_03.00/TS1015_i03.01_DESeq2_plotCounts_Design-group_GOI_2016.10.31") ## to change
#
getwd()
#
# notes
# Can run prev scripts or 
# Alternatively read in output from step a03.03 (saved on 2017.08.07)
load(file = 
"/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_R_Scripts_wd_03.00/TS1015_h03.04_DESeq2_DEA_Design-group_2016.10.31/Environment_TS1015_h03.03_DEA_setup_Design-group_2016.10.31.RData")
#
# Example code:
#_MEs= load(file = "./Environment_TS1015_a03.03_DEA_setup_tm.trt.int_2017.08.07.RData")
#


library("DESeq2") # start R session and load the DESeq2 package
library("ggplot2") # figures
library(dplyr) # for pipe function?
library(magrittr)
library(tidyr) # for gather function
library(knitr) # to combine? http://yihui.name/knitr/


#?plotCounts
colData(dds)


# Define the genes of interest
#--------------------------------to change
goi <- "Aqu2.1.32109_001, Aqu2.1.37571_001"
#-----------------------------------------


#Split a character vector that contains comma-separated values.
#data = '1.21, 1.985, 1.955, 2.015, 1.885';
# C = strsplit(data,', ')
# ref: http://stackoverflow.com/questions/8464312/convert-comma-separated-entry-to-columns
##################################
goi <- strsplit(goi, ', ') [[1]]
##################################

goi


# Join & tidy
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+0.5))) %>%
  	merge(colData(dds), ., by="row.names") %>%
  	gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

tcounts %>% 
  select(Row.names, time, treatment, gene, expression) %>% 
  head %>% 
  knitr::kable()

#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----
# single faceted plot  
pdf("TS1015_Cry2_ParA_splsda_2017.08.07-single-faceted-plot_forkeeps.pdf")  
for (i in goi) {
  p <- ggplot(filter(tcounts, gene==i), aes(time, expression, fill=treatment)) + 
  geom_boxplot() + 
  scale_fill_manual(values=col.list) +
  facet_wrap(~gene, scales="fixed", nrow=4, ncol=4) + 
  labs(x="Time (hours post release)", 
       y="Expression (log normalized counts)", 
       fill="Treatment", 
       title="TS1015_Cry2_ParA_splsda_2017.08.07")
  print(p)
}
dev.off()


pdf("TS1015_Cry2_ParA_splsda_2017.08.07-single-faceted-plot-2_forkeeps.pdf")       
ggplot(tcounts, aes(time, expression, fill=treatment)) + 
geom_boxplot() + 
scale_fill_manual(values=col.list) +
facet_wrap(~gene + treatment, scales="free_y") + 
labs(x="Time (hours post release)", 
	y="Expression (log normalized counts)", 
      fill="Treatment", 
      title="TS1015_Cry2_ParA_splsda_2017.08.07")
dev.off()
       
?facet_wrap  
# http://docs.ggplot2.org/0.9.3.1/facet_grid.html  
       


#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----
col.list <- c("gray98","dodgerblue4", "gray98","dodgerblue4", "gray98","dodgerblue4") # to change

tcounts$group <- factor(tcounts$group, levels = c("t1lgt", "t9lgt", "t1nat", "t9nat", "t1drk", "t9drk"))
# t1drk t1lgt t1nat t9drk t9lgt t9nat
# dds$treatment <- factor(dds$treatment, levels=c("nat","lgt", "drk"))
#x$name <- factor(x$name, levels = x$name[order(x$val)])
#x$name  # notice the changed order of factor levels  

pdf("TS1015_Cry2_ParA_splsda_2017.08.07-single-faceted-plot_forkeeps.pdf") 

ggplot(tcounts, aes(group, expression, fill=time), position=position_dodge(width=0)) + 
  geom_boxplot(outlier.size = 1, lwd=0.4, fatten = 1.5) + # lwd = line thickness, flatten = median; ref: http://stackoverflow.com/questions/23433776/change-thickness-of-the-whole-line-geom-boxplot.
  scale_fill_manual(values=col.list, breaks=c("lgt", "nat", "drk"), labels=c("Constant light", "Natural day-night cycle","Constant dark")) +
    # scale_fill_discrete(breaks=c("ctrl", "trt1", "trt2")labels=c("Control", "Treatment 1", "Treatment 2")) +
  #facet_wrap(~gene, scales="fixed") +
  facet_wrap(~gene, scales="free_y") +  
  labs(x="Developmental Time (hours post release)", 
       y="Gene Expression (log2 normalized counts)", 
       fill="Treatment" 
       #color="Percentile",
       #title="TS1015_Cry2_ParA_splsda_2017.08.07"
       ) +  
  guides(fill = guide_legend(order = 1), color = guide_legend(order = 2)) +

#scale_x_discrete(labels = c(quote(NA), quote(t1), quote(NA), quote(NA), quote(t9), quote(NA)), breaks =c("t1lgt", "t9lgt")) +
#scale_x_discrete(labels = c(quote(NA), quote(t1), quote(NA), quote(NA), quote(t9), quote(NA))) +
#########scale_x_discrete(labels = c("t1lgt", "t9lgt", "t1nat", "t9nat", "t1drk", "t9drk")) +
#####scale_y_continuous(limits = c(2, 13), breaks = c(0, 2, 4, 6, 8, 10, 12)) + # line spacing
#scale_y_continuous(limits =c(4, 10), breaks = c(4, 6, 8, 10)) + 
#scale_y_continuous(limits =c(2, 12), breaks = c(2, 4, 6, 8, 10, 12)) + 

# theme_classic() + # White background, no gridlines 
  theme_bw() + # White background with grid lines
  theme(   # remove the vertical grid lines
           # panel.grid.major.x =  element_line( size=.1, color="grey" ),
           # explicitly set the horizontal lines (or they will disappear too)
           panel.grid.major.y = element_line( size=.1, color="gray88" ),
           axis.text.y = element_text(size = 16), 
		axis.ticks = element_blank(), 
		panel.grid.major.x = element_blank(),
		#axis.text.x = element_text(size = 14, angle = 90),
		axis.text.x = element_blank(),
		axis.title.y = element_text(size = 16, angle = 90), #rel(0.7) ie x
		axis.title.x = element_text(size = 16),
		#title = element_text(size = 14),
           legend.justification =  c("right", "top")) 
           #legend.text = c(labels)) +         
# theme_minimal() # Minimal theme - no border
# theme_grey() # grey background (default theme)
# theme(axis.title.x = element_text(face="bold", colour="#990000", size=20),
#  axis.text.x  = element_text(angle=90, vjust=0.5, size=16)) +
  #geom_point(data=qt, aes(x=group, y=log2counts, color=log2counts), size = 0.5) + 
  #sc #+ scale_fill_discrete(breaks=c("0.25", "0.5", "0.75", "0.9", "0.95"), #labels=c("0.25", "0.5", "0.75", "0.9", "0.95"))

dev.off()


theme(axis.text.y = element_text(size = 3), 
		axis.ticks = element_blank(), 
		panel.grid.major.x = element_blank(),
		axis.text.x = element_text(size = 3, angle = 90),
		axis.title.y = element_text(size = 6, angle = 90), #rel(0.7) ie x
		axis.title.x = element_text(size = 6),
		title = element_text(size = 4),



#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----
pdf("TS1015_Cry2_ParA_splsda_2017.08.07-single-faceted-plot_forkeeps_quartiles.pdf") 

ggplot(tcounts, aes(group, expression, fill=treatment), position=position_dodge(width=0)) + 
  geom_boxplot(outlier.size = 1, lwd=0.4, fatten = 1.5) + # lwd = line thickness, flatten = median; ref: http://stackoverflow.com/questions/23433776/change-thickness-of-the-whole-line-geom-boxplot.
  scale_fill_manual(values=col.list, breaks=c("drk", "lgt", "nat"), labels=c("Dark", "Light", "Natural")) +
    # scale_fill_discrete(breaks=c("ctrl", "trt1", "trt2")labels=c("Control", "Treatment 1", "Treatment 2")) +
  facet_wrap(~gene, scales="fixed") + 
  labs(x="Developmental Time (hours post release)", 
       y="Gene Expression (log2 normalized counts)", 
       fill="Treatment", 
       color="Percentile",
       title="TS1015_Cry2_ParA_splsda_2017.08.07") +  
  guides(fill = guide_legend(order = 1), color = guide_legend(order = 2)) +

#scale_x_discrete(labels = c(quote(NA), quote(t1), quote(NA), quote(NA), quote(t9), quote(NA)), breaks =c("t1lgt", "t9lgt")) +
scale_x_discrete(labels = c(quote(NA), quote(t1), quote(NA), quote(NA), quote(t9), quote(NA))) +
scale_y_continuous(limits = c(2, 12), breaks = c(2, 3, 4, 5, 6, 7, 8, 10)) + # line spacing
#scale_y_continuous(limits =c(4, 10), breaks = c(4, 6, 8, 10)) + 
#scale_y_continuous(limits =c(2, 12), breaks = c(2, 4, 6, 8, 10, 12)) + 

# theme_classic() + # White background, no gridlines 
  theme_bw() + # White background with grid lines
  theme(   # remove the vertical grid lines
           panel.grid.major.x = element_blank(),
           
           # panel.grid.major.x =  element_line( size=.1, color="grey" ),

           # explicitly set the horizontal lines (or they will disappear too)
           panel.grid.major.y = element_line( size=.1, color="gray88" ),
           legend.justification =  c("right", "top")) +
           #legend.text = c(labels)) +         
# theme_minimal() # Minimal theme - no border
# theme_grey() # grey background (default theme)
# theme(axis.title.x = element_text(face="bold", colour="#990000", size=20),
#  axis.text.x  = element_text(angle=90, vjust=0.5, size=16)) +
  geom_point(data=qt, aes(x=group, y=log2counts, color=log2counts), size = 0.5) + 
  sc #+ scale_fill_discrete(breaks=c("0.25", "0.5", "0.75", "0.9", "0.95"), labels=c("0.25", "0.5", "0.75", "0.9", "0.95"))
dev.off()

# ----------------
# Figure legend    
# ----------------
#
#Log2 normalised gene expression levels in A. queenslandica larvae for <GENE NAME/ FAMILY NAME> which were identified as differentially expressed between two developmental time points (1 and 9 hours post release from the maternal sponge).  Alternating grey, yellow and green bars denote three different lighting treatments at each time point (constant dark, constant light and a natural day-night cycle, respectively).  Dashed lines show the transcriptome-wide percentiles (25th, 50th, 75th, 90th and 95th ) of transcript abundance in each treatment (after the initial filtering of lowly expressed genes). A doubling (or a reduction to 50%) is often considered as a biologically relevant change in gene expression, and this translates to one unit (+1 or -1) on the log2 scale (and this can be easily identified from the solid lines appearing on the x axis).
#


#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----
trawcounts <- t((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+0.5)) %>%
  	merge(colData(dds), ., by="row.names") %>%
  	gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

trawcounts %>% 
  select(Row.names, time, treatment, gene, expression) %>% 
  head %>% 
  knitr::kable()

# ----------------------------------
# single faceted plot  
pdf("TS1015_Cry2_ParA_splsda_2017.08.07-raw.counts-single-faceted-plot_forkeeps.pdf")  
  ggplot(trawcounts, aes(time, expression, fill=treatment)) + 
  geom_boxplot() + 
  scale_fill_manual(values=col.list) +
  facet_wrap(~gene, scales="fixed") + 
  labs(x="Time (hours post release)", 
       y="Expression (normalized counts)", 
       fill="Treatment", 
       title="TS1015_Cry2_ParA_splsda_2017.08.07")
dev.off()


# multi faceted plot  
pdf("TS1015_Cry2_ParA_splsda_2017.08.07-multi-ggplot2-catalog_normalised_counts.pdf")
for (i in goi) {
  p <- ggplot(filter(trawcounts, gene==i), aes(time, expression, color =factor(treatment))) + 
geom_boxplot() + 
  scale_fill_manual(values=col.list) +
  facet_wrap(~gene, scales="fixed") + 
  labs(x="Time (hours post release)", 
       y="Expression (normalized counts)", 
       fill="Treatment", 
       title="TS1015_Cry2_ParA_splsda_2017.08.07")
  ggtitle(i)
  print(p)
}
dev.off()

pdf("TS1015_Cry2_ParA_splsda_2017.08.07-raw.counts-single-faceted-plot-2_forkeeps.pdf")       
ggplot(trawcounts, aes(time, expression, fill=treatment)) + 
geom_boxplot() + 
scale_fill_manual(values=col.list) +
facet_wrap(~gene + treatment, scales="free_y") + 
labs(x="Time (hours post release)", 
	y="Expression (normalized counts)", 
      fill="Treatment", 
      title="TS1015_Cry2_ParA_splsda_2017.08.07")
dev.off()

#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----








