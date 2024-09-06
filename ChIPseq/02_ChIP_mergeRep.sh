#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 50G # The memory needed depends on the size of the genome and the size of fastqs
#SBATCH --cpus-per-task 1 # This allows to speed the mapping part of the pipeline
#SBATCH --time 04:00:00 # This depends on the size of the bedgraphs
#SBATCH --array=1-4 # Put here the rows from the table that need to be processed in the table
#SBATCH --job-name ChIPseq_mergeRep # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /scratch/ldelisle/Sheth/ChIPseq/ # This directory must exist, this is where will be the error and out files
## Specific to baobab:
##SBATCH --partition=shared-cpu # shared-cpu for CPU jobs that need to run up to 12h, public-cpu for CPU jobs that need to run between 12h and 4 days

gitHubDirectory=$1

# This script merge coverage of replicates.
# It uses deeptools bigwigAverage


##################################
#### TO SET FOR EACH ANALYSIS ####
##################################

binsize=50 # 50 is default value. Decreasing will increase computation time.

### Specify the paths (same as for script 01)

# Put in dirPathWithResults the directory
# where a directory will be created
# for each sample
dirPathWithResults="$PWD/"
# All merges are registered into a table where
# first column is the sample name
# second column is the name of all replicates separated by comma
filePathForTable="${gitHubDirectory}/ChIPseq/table_rep.txt"

### Specify the way to deal with dependencies:

# The only dependency is deeptools version >=3.5.4

# Or if you are using conda
# module purge
# module load Miniconda3/4.9.2

# You can choose to use a conda environment to solve deeptools>=3.5.4 dependency
# You can create it with: conda create -n deeptools3.5.5 deeptools=3.5.5
# Comment it if you will use module load
condaEnvName=deeptools3.5.5
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
python --version
if [ $? -ne 0 ]
then
  echo "python is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
bigwigAverage --version
if [ $? -ne 0 ]
then
  echo "deeptools is not installed or version is below 3.5.2. Please install it for example in the conda environment."
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

# Variable is initialized:
allBW=""

for s in $samples; do
    if [ -e ${dirPathWithResults}/${s}/${s}_macs_norm.bw ]; then
        allBW="${allBW} ${dirPathWithResults}/${s}/${s}_macs_norm.bw"
    fi
done

bigwigAverage --bigwigs $allBW --binSize $binsize -o ${sample}_norm.bw

mkdir -p ${dirPathWithResults}/allFinalFiles/bw_peaks
cp ${sample}_norm.bw ${dirPathWithResults}/allFinalFiles/bw_peaks/

mkdir -p ${dirPathWithResults}/toGEO
cp ${sample}_norm.bw ${dirPathWithResults}/toGEO/${sample}.bw
