#!/bin/bash

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 10G
#SBATCH --time 2:00:00
#SBATCH --job-name Rscripts

gitHubDirectory=$1
pathWithDependencies=$2


module purge
module load gcc/11.3.0
module load r/4.1.3


# This script will merge the tables 
path=/scratch/ldelisle/Sheth/RNAseq/

# time-course
configFile=${gitHubDirectory}/RNAseq/config/configFileRNAseq_step1_time_course.R
samplesPlan=${gitHubDirectory}/RNAseq/samplesPlan_time_course.txt
mkdir -p $(dirname $configFile)
echo "
### Required for all steps ###
RNAseqFunctionPath<-\"${pathWithDependencies}/rnaseq_rscripts/RNAseqFunctions.R\"
samplesPlan<-\"${samplesPlan}\" #This file should be a tabulated file with at least one column called sample. Optionnaly, the paths to the counts tables and FPKM tables can be provided under the column called: htseq_count_file and cufflinks_file.


#### STEP 1 - MERGE TABLES ###
#If the merged tables are not already generated:
outputFolderForStep1<-\"${path}/mergedTables/\"
#Needed for DESeq2:
mergeCounts<-T #Do you want to merge counts? T=yes F or commented=no
#Optional: subset the count table
subsetCounts<-T #Do you want to remove some genes from the count table
genesToRmFromCounts<-\"${path}/XYMTgenes.txt\" #List of genes id to remove (one per line with no header).
#Optional:
mergeFPKM<-T
oneLinePerEnsemblID<-T #By default cufflinks split the transcripts which do not overlap in different locus and so different lines, put T if you want to sum the FPKM for non overlapping transcripts (put F if not).
normFPKMWithAnoukMethod<-T #Anouk method: Genes that have the less variable rank should have the same expression.
chrToRemoveBeforeNormWithAnoukMethod<-c(\"chrX\",\"chrY\",\"chrM\")
" > ${configFile}

if [ ! -e ${path}/XYMTgenes.txt ]; then
  pathForGtf=$(ls ${path}/mer*UCSC.gtf)
  Rscript ${pathWithDependencies}toolBoxForMutantAndWTGenomes/scripts/getGeneListFromChrAndGTF.R $pathForGtf chrX,chrY,chrM ${path}
  mv ${path}/genesInchrX,chrY,chrMfrom* ${path}/XYMTgenes.txt
fi
# Adjust the samplesplan.txt
samplesPlanFile=$(cat $configFile| awk -F '=|<-|#' '$1=="samplesPlan"{print $2}' | tr -d "\"" | tr -d " ")
if [ ! $(grep "htseq_count_file" $samplesPlanFile) ]; then
  cat $samplesPlanFile | awk -v pa=$path 'BEGIN{print "htseq_count_file"}NR>1{print pa"/allFinalFiles/htseqCount_"$1".txt"}' > htseqCol.txt
  paste -d "\t" $samplesPlanFile htseqCol.txt > ${samplesPlanFile}_withPaths
  rm htseqCol.txt
fi
if [ ! $(grep "cufflinks_file" $samplesPlanFile) ]; then
  cat $samplesPlanFile | awk -v pa=$path 'BEGIN{print "cufflinks_file"}NR>1{print pa"/allFinalFiles/FPKM_"$1".txt"}' > CuffCol.txt
  if [ -e ${samplesPlanFile}_withPaths ];then
    mv ${samplesPlanFile}_withPaths tmp
    paste -d "\t" tmp CuffCol.txt > ${samplesPlanFile}_withPaths
    rm tmp
  else
    paste -d "\t" $samplesPlanFile CuffCol.txt > ${samplesPlanFile}_withPaths
  fi
  rm CuffCol.txt
fi
if [ -e ${samplesPlanFile}_withPaths ]; then
  cat $configFile | sed "s#$samplesPlanFile#${samplesPlanFile}_withPaths#" > ${configFile}_withPaths
  configFile=${configFile}_withPaths
fi
curWD="$PWD"
cd ${pathWithDependencies}rnaseq_rscripts/
echo "Version of rnaseq_rscripts"
git rev-parse HEAD
cd $curWD
# First generate tables
if [ ! -e ${path}/mergedTables/AllCufflinks_Simplified.txt ]; then
  Rscript ${pathWithDependencies}rnaseq_rscripts/step1-generateTables.R $configFile
  # copy to GEO:
  mkdir -p ${path}/toGEO
  cp ${path}/mergedTables/AllCufflinks_Simplified.txt ${path}/toGEO/AllCufflinks_Simplified_timecourse.txt
  cp ${path}/mergedTables/AllHTSeqCounts_subset.txt ${path}/toGEO/AllHTSeqCounts_subset_timecourse.txt
fi

if [ ! -e ${path}/DESeq2_pairwise/summary_significant.txt ]; then
  Rscript ${gitHubDirectory}/RNAseq/DESeq_pairwise_summary.R ${gitHubDirectory} ${path}
fi
