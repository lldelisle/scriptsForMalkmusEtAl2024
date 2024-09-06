
### Required for all steps ###
RNAseqFunctionPath<-"/home/ldelisle/softwares//rnaseq_rscripts/RNAseqFunctions.R"
samplesPlan<-"/home/ldelisle/softwares/scriptsForMalkmusEtAl2024//RNAseq/samplesPlan_time_course.txt" #This file should be a tabulated file with at least one column called sample. Optionnaly, the paths to the counts tables and FPKM tables can be provided under the column called: htseq_count_file and cufflinks_file.


#### STEP 1 - MERGE TABLES ###
#If the merged tables are not already generated:
outputFolderForStep1<-"/scratch/ldelisle/Sheth/RNAseq//mergedTables/"
#Needed for DESeq2:
mergeCounts<-T #Do you want to merge counts? T=yes F or commented=no
#Optional: subset the count table
subsetCounts<-T #Do you want to remove some genes from the count table
genesToRmFromCounts<-"/scratch/ldelisle/Sheth/RNAseq//XYMTgenes.txt" #List of genes id to remove (one per line with no header).
#Optional:
mergeFPKM<-T
oneLinePerEnsemblID<-T #By default cufflinks split the transcripts which do not overlap in different locus and so different lines, put T if you want to sum the FPKM for non overlapping transcripts (put F if not).
normFPKMWithAnoukMethod<-T #Anouk method: Genes that have the less variable rank should have the same expression.
chrToRemoveBeforeNormWithAnoukMethod<-c("chrX","chrY","chrM")

