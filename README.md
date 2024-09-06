# scriptsForMalkmusEtAl2024

All scripts necessary to build figures from raw data in Malkmus et al. 2024.

Each of these directories contain command lines for the corresponding analyses: `ChIPseq`, `RNAseq`. They also contains figures with quantifications in the `outputs` folder.

The `plot` directory contains command line for the [pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks) plots with ini files and output pdf files.

## Run all analyzes from scratch

### What has been done before

The fasta file for mm10 has been downloaded. The STAR index has been generated before we start but it needs to be at the good place.

```bash
cp -r /work/updub/scratch/ldelisle/genomes/STARIndex_2.7.0e/* /scratch/ldelisle/genomes/STARIndex_2.7.0e/
```

The bowtie2 index has been generated.

The conda environments have been created
```bash
conda create -n rnaseq_2021 python=3.6 mamba
conda activate rnaseq_2021
mamba install cutadapt=1.16 star=2.7.0e samtools=1.9 bedtools=2.31.1 ucsc-bedgraphtobigwig=357
conda deactivate

conda create -n deeptools3.5.5 deeptools=3.5.5

conda create -n atac202209 cutadapt=4.1 bowtie2=2.4.5 bedtools=2.30.0 samtools=1.16.1 macs2=2.2.7.1 ucsc-bedgraphtobigwig=377 picard=2.27.4
```

### RNAseq

The directory is created and fastqs are uploaded:

```bash
mkdir -p /scratch/ldelisle/Sheth/RNAseq/fastq/
```

Run the first step:
```bash
sbatch RNAseq/01_RNAseq.sh $PWD/
```

Run the second:
```bash
sbatch --dependency=afterok:24057022 RNAseq/02_RNAseq_mergeRep.sh $PWD/
```

Clone the 2 repos: https://github.com/lldelisle/rnaseq_rscripts and https://github.com/lldelisle/toolBoxForMutantAndWTGenomes in a common directory (in my case ~/softwares/).

Prepare the directory for config files

```bash
mkdir -p RNAseq/config/
```

Run the third sbatch script

```bash
sbatch --dependency=afterok:24057022 --chdir RNAseq/ RNAseq/03_run_global_Rscripts.sh $PWD/ ~/softwares/
```

Run locally make_figure.R
```bash
Rscript RNAseq/make_figures.R
```

### ChIPseq

The directory is created and fastqs are uploaded:

```bash
mkdir -p /scratch/ldelisle/Sheth/ChIPseq/fastq/
```

The first script is run:
```bash
sbatch ChIPseq/01_ChIP_PE.sh $PWD/
```

The second script is run:
```bash
sbatch --dependency=afterok:24061502 ChIPseq/02_ChIP_mergeRep.sh $PWD/
```

The third script is run:
```bash
sbatch --dependency=afterok:24061502 ChIPseq/03_ChIP_PE_consensus_peaks.sh $PWD/
```

Some consensus peaks are copied back to the git repo:
```bash
mkdir -p ChIPseq/outputs/
cp /scratch/ldelisle/Sheth/ChIPseq/toGEO/b-Catenin_*DMSO*consensus.narrowPeak ChIPseq/outputs/
gzip ChIPseq/outputs/b-Catenin_*DMSO*consensus.narrowPeak
```

Locally I run the [fourth bash script](./ChIPseq/04_deeper_analysis.sh) line by line.

Run the fifth script
```bash
sbatch ChIPseq/05_deeptools.sh $PWD/
```

Run HOMER analysis on Galaxy, the command line and explanations are in `ChIPseq/06_homer_motif.sh`.

### pyGenomeTracks

Run locally
```bash
bash plots/plots.sh
```
