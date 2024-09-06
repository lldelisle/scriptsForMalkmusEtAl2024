#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 50G # The memory needed depends on the size of the genome and the size of fastqs
#SBATCH --cpus-per-task 1 # This allows to speed the mapping part of the pipeline
#SBATCH --time 04:00:00 # This depends on the size of the bedgraphs
#SBATCH --array=1-4 # Put here the rows from the table that need to be processed in the table
#SBATCH --job-name ChIPseq_consensusPeaks # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /scratch/ldelisle/Sheth/ChIPseq/ # This directory must exist, this is where will be the error and out files
## Specific to baobab:
##SBATCH --partition=shared-cpu # shared-cpu for CPU jobs that need to run up to 12h, public-cpu for CPU jobs that need to run between 12h and 4 days

gitHubDirectory=$1

# This script generate a narrowPeak consensus between replicates


##################################
#### TO SET FOR EACH ANALYSIS ####
##################################
# Number of replicates where a peak should be present
minNumberOverlap=2
# Genome size
# Used by MACS2: H. sapiens: 2700000000, M. musculus: 1870000000, D. melanogaster: 120000000, C. elegans: 90000000
effectiveGSize=1870000000

### Specify the paths (same as for script 01)

# Put in dirPathWithResults the directory
# where a directory will be created
# for each sample
dirPathWithResults="$PWD/"
# All merges are registered into a table where
# first column is the sample name
# second column is the name of all replicates separated by comma
filePathForTable="${gitHubDirectory}/ChIPseq/table_rep.txt"

picardCommand="picard"
# For baobab:
# picardCommand="java -jar $EBROOTPICARD/picard.jar"

### Specify the way to deal with dependencies:
# Dependencies are samtools, macs2, picard, bedtools

# Or if you are using conda
# module purge
# module load Miniconda3/4.9.2

# Comment it if you will use module load
condaEnvName=atac202209
######


##################################
####### BEGINING OF SCRIPT #######
##################################

# Check everything is set correctly:
if [ ! -z ${condaEnvName} ]; then
    # Conda environment:
    # This line is to adapt the conda to the shell
    source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
    # We check if the conda environment exists
    exists=$(conda info --envs | awk -v ce=${condaEnvName} '$1==ce{print}' | wc -l)
    # It if does not exists an error is raised
    if [ $exists -ne 1 ]; then
    echo "conda environment ${condaEnvName} does not exists. Create it before."
    exit 1
    fi
    # Activate the conda environment
    conda activate ${condaEnvName}
fi

# Check all softwares are present and write version to stdout:
samtools --version
if [ $? -ne 0 ]
then
  echo "samtools is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
macs2 --version
if [ $? -ne 0 ]
then
  echo "macs2 is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
v=$($picardCommand MarkDuplicates --version 2>&1)
if [ $? -eq 127 ]
then
  echo "picard is not installed but required. Or picardCommand is wrongly set. Please install it for example in the conda environment."
  exit 1
fi
echo "picard MarkDuplicates version $v"
bedtools --version
if [ $? -ne 0 ]
then
  echo "bedtools is not installed but required. Please install it for example in the conda environment."
  exit 1
fi

sample=$(cat ${filePathForTable} | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}')
samples=$(cat ${filePathForTable} | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{split($2,a,",");for(j in a){print a[j]}}')
n=$(echo $samples | wc -w)

# Each sample is processed into an independent directory:
pathResults=${dirPathWithResults}/${sample}/

# The directory is created (if not existing)
mkdir -p ${pathResults}

# The name of the sample is written in stdout
echo $sample

# The analysis part takes part within the pathResults
cd $pathResults

# Generate a narrowPeak which uses as many reads of each replicates
# Get the number of reads per bam
nreads=""
for s in $samples; do
  if [ ! -e ${dirPathWithResults}/${s}/${s}_mapped_sorted_q30_rmdup.bam ]; then
    if [ ! -e ${dirPathWithResults}/${s}/${s}_mapped_sorted_q30.bam ]; then
      echo "cannot find bam for $s"
      exit 1
    fi
    ${picardCommand} MarkDuplicates SORTING_COLLECTION_SIZE_RATIO=0.15 I=${dirPathWithResults}/${s}/${s}_mapped_sorted_q30.bam O=${dirPathWithResults}/${s}/${s}_mapped_sorted_q30_rmdup.bam M=${dirPathWithResults}/${s}/${s}_mapped_sorted_q30_rmdup.log REMOVE_DUPLICATES=true AS=true
  fi
  nreads="$nreads $(samtools view -c ${dirPathWithResults}/${s}/${s}_mapped_sorted_q30_rmdup.bam)"
done

finalreads=$(echo $nreads | awk 'BEGIN{minval=1e15}{for(i=1;i<=NF;i++){if($i < minval){minval=$i}}}END{print minval}')

# Downsample bams
i=1
bams=""
for s in $samples; do
  if [ ! -e ${s}_mapped_sorted_q30_rmdup_sub.bam ]; then
    frac=$(echo $nreads | awk -v i=$i -v fr=${finalreads} '{print fr/$i}')
    if [ "$frac" = "1" ]; then
      cp ${dirPathWithResults}/${s}/${s}_mapped_sorted_q30_rmdup.bam ${s}_mapped_sorted_q30_rmdup_sub.bam
    else
      samtools view -b -o ${s}_mapped_sorted_q30_rmdup_sub.bam --subsample $frac --subsample-seed 1 ${dirPathWithResults}/${s}/${s}_mapped_sorted_q30_rmdup.bam
    fi
  fi
  bams="$bams ${s}_mapped_sorted_q30_rmdup_sub.bam"
  i=$((i + 1))
done

# Call peaks with all BAMs
if [ ! -e ${sample}_macs_peaks.narrowPeak ]; then
  macs2 callpeak -g $effectiveGSize -t $bams -n ${sample}_macs --call-summits --format BAMPE 2> ${sample}_macs.log
fi

# Compare this narrowPeak individual narrowPeaks
indivnp=""
for s in $samples; do
  if [ ! -e ${dirPathWithResults}/${s}/${s}_macs_peaks.narrowPeak ]; then
    echo "cannot find narrowPeak for $s"
    exit 1
  fi
  indivnp="$indivnp ${dirPathWithResults}/${s}/${s}_macs_peaks.narrowPeak"
done

# Intersect them
if [ ! -e multi_inter.bed ]; then
  bedtools multiinter -i $indivnp > multi_inter.bed
fi

# Filter
if [ ! -e multi_inter_filtered.bed ]; then
  cat multi_inter.bed | awk -v minNumberOverlap=$minNumberOverlap '$4>=minNumberOverlap{print}' > multi_inter_filtered.bed
fi

# Overlap
if [ ! -e ${sample}_macs_peaks_inter.txt ]; then
  bedtools intersect -a ${sample}_macs_peaks.narrowPeak -b multi_inter_filtered.bed \
    -wa -wb > ${sample}_macs_peaks_inter.txt
fi

# Select narrowPeaks with summit overlapping:
if [ ! -e ${sample}_consensus.narrowPeak ]; then
  cat ${sample}_macs_peaks_inter.txt | awk '($2+$10)>=$12 && ($2+$10)<$13{print}' | cut -f 1-10 | uniq > ${sample}_consensus.narrowPeak
fi

mkdir -p ${dirPathWithResults}/allFinalFiles/bw_peaks
cp ${sample}_consensus.narrowPeak ${dirPathWithResults}/allFinalFiles/bw_peaks/

mkdir -p ${dirPathWithResults}/toGEO
cp ${sample}_consensus.narrowPeak ${dirPathWithResults}/toGEO/
