#!/bin/bash

#SBATCH -o slurm-%x-%A.out # Template for the std output of the job uses the job name, the job id
#SBATCH -e slurm-%x-%A.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 10G # The memory needed depends on the size of the genome and the size of fastqs
#SBATCH --cpus-per-task 4 # This allows to speed the pipeline
#SBATCH --time 04:00:00 # This depends on the size of the bedgraphs
#SBATCH --job-name deeptools # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /scratch/ldelisle/Sheth/ChIPseq/ # This directory must exist, this is where will be the error and out files
## Specific to baobab:
##SBATCH --partition=shared-cpu # shared-cpu for CPU jobs that need to run up to 12h, public-cpu for CPU jobs that need to run between 12h and 4 days

gitHubDirectory=$1

dirPathWithResults=$PWD/
binsize=1000
condaEnvName=deeptools3.5.5
sample=FL_E10.5

module purge
module load gcc/11.3.0 samtools/1.14

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

mkdir -p deeptools
cd deeptools

# Plot the heatmap
bws=""
labels=""
for prot in b-Catenin H3K27Ac; do
    for treatment in DMSO C59; do
        fullsample=${prot}_${sample}+${treatment}+6h
        bws="$bws ${dirPathWithResults}/${fullsample}/${fullsample}_norm.bw"
        labels="$labels $fullsample"
    done
done
if [ ! -e bcat_${sample}.npz ]; then
    computeMatrix reference-point \
        --regionsFileName \
            ${gitHubDirectory}/ChIPseq/outputs/bcat_with_TSS_${sample}.bed \
            ${gitHubDirectory}/ChIPseq/outputs/bcat_noTSS_${sample}.bed \
        --scoreFileName ${bws} --outFileName bcat_${sample}.npz \
        --samplesLabel ${labels} --numberOfProcessors 4 \
        --outFileNameMatrix bcat_${sample}.txt \
        --beforeRegionStartLength 2500 --afterRegionStartLength 2500 \
        --sortRegions 'descend' --sortUsing 'mean' \
        --averageTypeBins 'mean' --missingDataAsZero \
        --binSize 10 --referencePoint center
fi
plotHeatmap --matrixFile bcat_${sample}.npz \
    --outFileName ${gitHubDirectory}/ChIPseq/outputs/heatmap_${sample}.pdf  \
    --plotFileFormat 'pdf' --outFileSortedRegions bcat_${sample}_regions.txt \
    --dpi '200' --sortRegions 'descend' --sortUsing 'mean' \
    --averageTypeSummaryPlot 'mean' --plotType 'lines' \
    --missingDataColor 'black'   --alpha '1.0' \
    --sortUsingSamples 1 --xAxisLabel 'distance from bCat summit' \
    --yAxisLabel 'coverage'  --heatmapWidth 7.5 --heatmapHeight 25.0 \
    --whatToShow 'plot, heatmap and colorbar' --refPointLabel 'summit' \
    --plotTitle "Heatmap for bCat ${sample} DMSO (more or less than 1kb from TSS)" \
    --legendLocation 'best'  --labelRotation '0'

# Generate the correlation matrix
bams=""
labels=""
for prot in b-Catenin H3K27Ac; do
    for treatment in DMSO C59; do
        for rep in 1 2; do
            fullsample=${prot}_${sample}+${treatment}+6h_rep${rep}
            if [ ! -e ${dirPathWithResults}/${fullsample}/${fullsample}_mapped_sorted_q30_rmdup.bam ]; then
                echo "cannot find rmdup bam for $fullsample"
                exit 1
            fi
            if [ ! -e ${dirPathWithResults}/${fullsample}/${fullsample}_mapped_sorted_q30_rmdup.bam.bai ]; then
                samtools index ${dirPathWithResults}/${fullsample}/${fullsample}_mapped_sorted_q30_rmdup.bam
            fi
            bams="$bams ${dirPathWithResults}/${fullsample}/${fullsample}_mapped_sorted_q30_rmdup.bam"
            labels="$labels $fullsample"
        done
    done
done
if [ ! -e counts_${binsize}_${sample}.npz ]; then
    multiBamSummary bins \
        --outFileName counts_${binsize}_${sample}.npz --bamfiles $bams \
        --labels $labels --outRawCounts raw_counts_${binsize}_${sample}.txt \
        --binSize ${binsize} --distanceBetweenBins '0' \
        --numberOfProcessors 4
fi
# Plot the correlation
plotCorrelation --corData counts_${binsize}_${sample}.npz --plotFile ${gitHubDirectory}/ChIPseq/outputs/correlation_${binsize}_${sample}.pdf --corMethod 'pearson' --whatToPlot 'heatmap' --colorMap 'RdYlBu_r' --plotNumbers --plotTitle 'Correlation on $((binsize / 1000))kb bins genome-wide'  --plotWidth 11.0 --plotHeight 9.5 --skipZeros --plotFileFormat 'pdf' --removeOutliers --outFileCorMatrix ${gitHubDirectory}/ChIPseq/outputs/correlation_${binsize}_${sample}.txt
