#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 50G
#SBATCH --cpus-per-task 16
#SBATCH --time 10:00:00
#SBATCH --array=1-30
#SBATCH --job-name RNAseq_PE
#SBATCH --chdir /scratch/ldelisle/Sheth/RNAseq/

gitHubDirectory=$1

# This script run cutadapt to remove adapters and bad quality bases
# make alignment with STAR ENCODE parameters
# Evaluate FPKM with cufflinks
# Coverage normalized to million mapped reads with STAR

path="$PWD/"
pathForFastq="$path/fastq/"
pathForTable="${gitHubDirectory}/RNAseq/table.txt"
pathForFasta="/home/ldelisle/genomes/fasta/"
pathForIndex="/scratch/ldelisle/genomes/STARIndex_2.7.0e/"
# To copy from /work/updub/scratch/ldelisle/genomes/STARIndex_2.7.0e/
genome=mm10
ensemblVersion=102
versionOfCufflinks="2.2.1"
# I tryed 36 and it failed:
# BAMoutput.cpp:27:BAMoutput: exiting because of *OUTPUT FILE* error: could not create output file ./_STARtmp//BAMsort/19/48
# SOLUTION: check that the path exists and you have write permission for this file. Also check ulimit -n and increase it to allow more open files.
nbOfThreads=16

condaEnvName=rnaseq_2021

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

# cutadapt version 1.16
# star version 2.7.0e
# samtools version 1.9
# bedtools 2.31.1
# ucsc-bedgraphtobigwig 357

sample=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}')
fastqR1File=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')
fastqR2File=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $3}')
stranded=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $4}')
adapterSeq="TruSeq"

if [ ! -z $newGenome ]; then
  genome=$newGenome
fi

if [ $genome = "mm10" ]; then
  versionOfGtf="Mus_musculus.GRCm38.$ensemblVersion"
  gtfFile=${path}mergeOverlapGenesOfFilteredTranscriptsOf${versionOfGtf}_ExonsOnly_UCSC.gtf
else
  echo "unknown genome"
  exit 1
fi
indexPath=${pathForIndex}${genome}

#For the first one, The gtf is downloaded.
if [ $SLURM_ARRAY_TASK_ID == 1 ] && [ ! -e ${pathForFasta}${genome}.fa ];then
  cat ${pathForFasta}$genome/*.fa.gz > ${pathForFasta}${genome}.fa.gz
  gunzip ${pathForFasta}${genome}.fa.gz
fi
if [ $SLURM_ARRAY_TASK_ID == 1 ] && [ ! -e $gtfFile ];then
  wget "https://zenodo.org/record/4596490/files/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsOnly_UCSC.gtf.gz?download=1" -O ${gtfFile}.gz
  gunzip ${gtfFile}.gz
fi
if [ $SLURM_ARRAY_TASK_ID == 1 ] && [ ! -e ${path}/MTmouse.gtf ];then
  echo -e "chrM\tchrM_gene\texon\t0\t16299\t.\t+\t.\tgene_id \"chrM_gene_plus\"; transcript_id \"chrM_tx_plus\"; exon_id \"chrM_ex_plus\";">MTmouse.gtf
  echo -e "chrM\tchrM_gene\texon\t0\t16299\t.\t-\t.\tgene_id \"chrM_gene_minus\"; transcript_id \"chrM_tx_minus\"; exon_id \"chrM_ex_minus\";" >>MTmouse.gtf
fi

mkdir -p ${path}STAR/$sample

pathResults=${path}STAR/${sample}/
echo $sample
cd $pathResults

if [ ! -e ${path}allFinalFiles/reports/${sample}_report-cutadapt_PE.txt ]; then
  fastqR1="${pathForFastq}/${fastqR1File}"
  fastqR2="${pathForFastq}/${fastqR2File}"
  if [ ! -e $fastqR1 ]; then
    mkdir -p $pathForFastq
    cd $pathForFastq
    module load sra-toolkit/2.9.6
    fasterq-dump -o ${sample}.fastq ${sra}
    mv ${sample}_1.fastq ${sample}_R1.fastq
    mv ${sample}_2.fastq ${sample}_R2.fastq
    gzip ${sample}_R1.fastq
    gzip ${sample}_R2.fastq
    cd $pathResults
    fastqR1="${pathForFastq}/${sample}_R1.fastq.gz"
    fastqR2="${pathForFastq}/${sample}_R2.fastq.gz"
  fi
  if [ $adapterSeq = "TruSeq" ]; then
    cutadapt -j $nbOfThreads -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -q 30 -m 15 -o ${pathResults}${sample}-cutadapt_R1.fastq.gz -p ${pathResults}${sample}-cutadapt_R2.fastq.gz $fastqR1 $fastqR2 > ${pathResults}${sample}_report-cutadapt_PE.txt
  else
    if [ $adapterSeq = "Nextera" ]; then
      cutadapt -j $nbOfThreads -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -q 30 -m 15 -o ${pathResults}${sample}-cutadapt_R1.fastq.gz -p ${pathResults}${sample}-cutadapt_R2.fastq.gz $fastqR1 $fastqR2 > ${pathResults}${sample}_report-cutadapt_PE.txt
    else
      echo "YOU NEED TO WRITE THE CODE"
      exit 1
    fi
  fi
  mkdir -p ${path}allFinalFiles/reports/
  cp ${pathResults}${sample}_report-cutadapt_PE.txt ${path}allFinalFiles/reports/
fi

if [ ! -e ${path}allFinalFiles/bam/${sample}_Aligned.sortedByCoord.out.bam ];then
  if [ "$stranded" = "unstranded" ]; then 
    #I need to add --outSAMstrandField intronMotif because it is not stranded library
    STAR --runThreadN $nbOfThreads --genomeDir ${indexPath} --readFilesIn ${pathResults}${sample}-cutadapt_R1.fastq.gz ${pathResults}${sample}-cutadapt_R2.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate  --outSAMstrandField intronMotif  --sjdbOverhang '99' --sjdbGTFfile $gtfFile  --quantMode GeneCounts  --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outWigType bedGraph --outWigStrand Unstranded
  else
    STAR --runThreadN $nbOfThreads --genomeDir ${indexPath} --readFilesIn ${pathResults}${sample}-cutadapt_R1.fastq.gz ${pathResults}${sample}-cutadapt_R2.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate  --sjdbOverhang '99' --sjdbGTFfile $gtfFile  --quantMode GeneCounts  --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outWigType bedGraph --outWigStrand Stranded
  fi
  mkdir -p ${path}allFinalFiles/reports/
  cp Log.final.out ${path}allFinalFiles/reports/${sample}_STAR_logFinal.txt
  mkdir -p ${path}allFinalFiles/bam/
  cp ${pathResults}Aligned.sortedByCoord.out.bam ${path}allFinalFiles/bam/${sample}_Aligned.sortedByCoord.out.bam
  samtools index ${path}allFinalFiles/bam/${sample}_Aligned.sortedByCoord.out.bam
fi

if [ -e ${pathForFasta}${genome}.fa ] && [ -e ${path}MTmouse.gtf ] && [ -e $gtfFile ] && [ -s ${pathResults}Aligned.sortedByCoord.out.bam ];then
  if [ ! -e ${path}allFinalFiles/FPKM_${sample}_isoforms.txt ];then
    echo "export PATH=$PATH:/home/ldelisle/softwares/cufflinks-${versionOfCufflinks}.Linux_x86_64" >cufflinks_${sample}.sh
    echo "mkdir -p ${pathResults}cufflinksWOMT" >>cufflinks_${sample}.sh
    echo "cufflinks -p 10 -o ${pathResults}cufflinksWOMT --max-bundle-length 10000000 --multi-read-correct --library-type \"fr-firststrand\" -b ${pathForFasta}${genome}.fa  --no-effective-length-correction -M ${path}MTmouse.gtf -G $gtfFile ${pathResults}Aligned.sortedByCoord.out.bam" >>cufflinks_${sample}.sh
    echo "" >>cufflinks_${sample}.sh
    echo "mkdir -p ${path}allFinalFiles" >>cufflinks_${sample}.sh
    echo "cp ${pathResults}cufflinksWOMT/genes.fpkm_tracking ${path}allFinalFiles/FPKM_${sample}.txt" >>cufflinks_${sample}.sh
    echo "cp ${pathResults}cufflinksWOMT/isoforms.fpkm_tracking ${path}allFinalFiles/FPKM_${sample}_isoforms.txt" >>cufflinks_${sample}.sh
    
    if [ "$stranded" = "unstranded" ]; then 
      sed -i 's/fr-firststrand/fr-unstranded/g' cufflinks_${sample}.sh
    fi
    echo "Launching cufflinks"
    bash cufflinks_${sample}.sh &
  fi
else
 echo "cufflinks not launch because some files are missing."
fi


# if { [ ! -e accepted_hits_unique_${sample}.bam ] && [ -s ${pathResults}Aligned.sortedByCoord.out.bam ] ;} || [ -e tmp.header ] ;then
#  echo "Compute uniquely aligned"
#  samtools view -H Aligned.sortedByCoord.out.bam >  tmp.header
#  samtools view -@ 5 Aligned.sortedByCoord.out.bam | grep  -w "NH:i:1" | cat tmp.header - | samtools view -@ 5 -b > accepted_hits_unique_${sample}.bam
#  rm tmp.header
# fi

if [ ! -e ${path}allFinalFiles/htseqCount_${sample}.txt ] && [ -e ReadsPerGene.out.tab ];then 
  mkdir -p ${path}allFinalFiles
  echo "write htseqCount"
  if [ "$stranded" = "unstranded" ]; then 
    cat ReadsPerGene.out.tab | awk '{print $1"\t"$2}' > ${path}allFinalFiles/htseqCount_${sample}.txt
  else
    cat ReadsPerGene.out.tab | awk '{print $1"\t"$4}' > ${path}allFinalFiles/htseqCount_${sample}.txt
  fi
fi


# Coverage
# Sort files
for f in *.out.bg; do
  output=${f/.bg/.sorted.bg}
  if [ ! -e $output ]; then
    bedtools sort -i $f > $output
  fi
done

# Convert files
if [ "$stranded" = "unstranded" ]; then 
  if [ ! -e ${sample}_Uniq_norm.bw ]; then
    bedGraphToBigWig Signal.Unique.str1.out.sorted.bg ${pathForFasta}${genome}.fa.fai ${sample}_Uniq_norm.bw
  fi
elif [ "$stranded" = "forward" ]; then
  if [ ! -e ${sample}_Uniq_fwd_norm.bw ]; then
    bedGraphToBigWig Signal.Unique.str1.out.sorted.bg ${pathForFasta}${genome}.fa.fai ${sample}_Uniq_fwd_norm.bw
  fi
  if [ ! -e ${sample}_Uniq_rev_norm.bw ]; then
  bedGraphToBigWig Signal.Unique.str2.out.sorted.bg ${pathForFasta}${genome}.fa.fai ${sample}_Uniq_rev_norm.bw
  fi
else
  if [ ! -e ${sample}_Uniq_fwd_norm.bw ]; then
    bedGraphToBigWig Signal.Unique.str2.out.sorted.bg ${pathForFasta}${genome}.fa.fai ${sample}_Uniq_fwd_norm.bw
  fi
  if [ ! -e ${sample}_Uniq_rev_norm.bw ]; then
  bedGraphToBigWig Signal.Unique.str1.out.sorted.bg ${pathForFasta}${genome}.fa.fai ${sample}_Uniq_rev_norm.bw
  fi
fi

mkdir -p ${path}/allFinalFiles/bw
cp *bw ${path}/allFinalFiles/bw/


wait
echo "Everything is done"
find . -size 0 -delete

mkdir -p ${path}/toGEO
cp *bw ${path}/toGEO/
