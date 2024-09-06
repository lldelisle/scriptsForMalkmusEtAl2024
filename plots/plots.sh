# Using pygenometracks version 3.9
gitHubDirectory=$PWD/
GEODirectory=~/mnt/scratch/Sheth/

condaEnvName=lastVersion # 3.9
# This line is to adapt the conda to the shell
source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
# Activate the conda environment
conda activate ${condaEnvName}

pathWithChIP="${GEODirectory}/ChIPseq/toGEO"
pathWithRNAseq="${GEODirectory}/RNAseq/toGEO"

ymax_values=('1' '3' '3' '4')
genes=('Fgf8' 'Axin2' 'Grem1' 'Bmp4')
scalebars=('20' '50' '10' '50')

mkdir -p ${gitHubDirectory}/plots/outputs
cd ${gitHubDirectory}/plots/outputs

for i in "${!ymax_values[@]}"; do
    ini_file=${genes[$i]}_Main.ini
    echo "[genes]
file = ${gitHubDirectory}/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz
color = black
color_utr = black
prefered_name = gene_name
merge_transcripts = true
merge_overlapping_exons = true
all_labels_inside = true
labels_in_margin = true
fontstyle = oblique
fontsize = 14
height = 3

[spacer]
height = 0.04

[CRMs]
file = ${gitHubDirectory}/plots/annotations.bed
display = collapsed
labels = false
color = red

[spacer]
height = 0.04

[bcat]
file = ${pathWithChIP}/b-Catenin_FL_E10.5+DMSO+6h.bw
color = black
height = 3
title = b-Catenin DMSO
nans_to_zeros = true
min_value = 0
max_value = ${ymax_values[$i]}

[bcat_C59]
file = ${pathWithChIP}/b-Catenin_FL_E10.5+C59+6h.bw
color = orange
height = 3
title = b-Catenin C59
nans_to_zeros = true
min_value = 0
max_value = ${ymax_values[$i]}
orientation = inverted

[scalebar]
file_type = scalebar
where = bottom
size = ${scalebars[$i]}000
" > ${ini_file}
done

pyGenomeTracks --tracks Axin2_Main.ini --region chr11:108350000-109350000 -o Axin2_Main.pdf

pyGenomeTracks --tracks Grem1_Main.ini --region chr2:113400000-113830000 -o Grem1_Main.pdf

pyGenomeTracks --tracks Bmp4_Main.ini --region chr14:45700000-46900000 -o Bmp4_Main.pdf

pyGenomeTracks --tracks Fgf8_Main.ini --region chr19:45425000-46000000 -o Fgf8_Main.pdf

# For the supplement, we need the CC from Andrey et al. 2017 lifted over on mm10
if [ ! -e CC_liftover ]; then
    tmpdir=$(mktemp -d)

    wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2473963&format=file&file=GSM2473963%5FCC%2DFL%2DE105%2DWt%2DMm%5FMerged%2DSmoothed%2D3kb%2DNorm%2Etar%2Egz" -O ${tmpdir}/CC.tar.gz

    tar zxvmf ${tmpdir}/CC.tar.gz -C ${tmpdir}

    # For liftover
    # Get liftOver for linux
    wget "https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver" -O "${tmpdir}/liftOver" -nc
    chmod +x ${tmpdir}/liftOver
    # Get the chainFile
    wget "https://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz" -O ${tmpdir}/mm9ToMm10.over.chain.gz -nc
    mkdir -p CC_liftover
    for gene in Axin2 Bmp4 Fgf8 Grem1 Shh; do
        awk '$0~/^chr/{print}' ${tmpdir}/CC-FL-E105-Wt-Mm_Merged-Smoothed-3kb-Norm/*${gene}.3kb.bedgraph > ${tmpdir}/${gene}.3kb.bedgraph
        ${tmpdir}/liftOver ${tmpdir}/${gene}.3kb.bedgraph ${tmpdir}/mm9ToMm10.over.chain.gz CC_liftover/${gene}.3kb.bedgraph ${tmpdir}/unmapped.${gene}.3kb.bedgraph
    done
fi

ymax_values_CC=('30' '45' '30' '40' '40')
ymax_values_RNA=('4' '3' '5' '8' '8')
ymax_values_bcat=('3' '3' '3' '1' '1')
ymax_values_k27=('4' '1' '3' '5' '5')
genes=('Axin2' 'Grem1' 'Bmp4' 'Fgf8' 'Shh')

for i in "${!genes[@]}"; do
    ini_file=${genes[$i]}_Supplement.ini
    echo "[genes]
file = ${gitHubDirectory}/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz
color = black
color_utr = black
prefered_name = gene_name
merge_transcripts = true
merge_overlapping_exons = true
all_labels_inside = true
labels_in_margin = true
fontstyle = oblique
fontsize = 14
height = 3

[spacer]
height = 0.04

[CRMs]
file = ${gitHubDirectory}/plots/annotations.bed
display = collapsed
labels = false
color = red

[CC]
file = CC_liftover/${genes[$i]}.3kb.bedgraph
title = CC E10.5 Andrey et al. 2017
min_value = 0
max_value = ${ymax_values_CC[$i]}
color = grey
height = 3 
use_middle = true

[spacer]
height = 0.04
" > ${ini_file}
for time in 6 12 18; do
    echo "[RNA_DMSO_${time}]
file = ${pathWithRNAseq}/RNA_FL_timecourse_E10.5+DMSO+${time}h_Uniq_norm.bw
title = RNA-seq Ctrl vs C59 E10.5+${time}h
nans_to_zeros = true
min_value = 0
max_value = ${ymax_values_RNA[$i]}
color = black
height = 3

[RNA_C59]
file = ${pathWithRNAseq}/RNA_FL_timecourse_E10.5+C59+${time}h_Uniq_norm.bw
nans_to_zeros = true
color = orange
show_data_range = false
overlay_previous = share-y

[spacer]
height = 0.04
" >> ${ini_file}
done
echo "[bcat]
file = ${pathWithChIP}/b-Catenin_FL_E10.5+DMSO+6h.bw
color = black
height = 3
title = b-Catenin DMSO
nans_to_zeros = true
min_value = 0
max_value = ${ymax_values_bcat[$i]}

[bcat_C59]
file = ${pathWithChIP}/b-Catenin_FL_E10.5+C59+6h.bw
color = orange
height = 3
title = b-Catenin C59
nans_to_zeros = true
min_value = 0
max_value = ${ymax_values_bcat[$i]}
orientation = inverted

[spacer]
height = 0.04

[ac]
file = ${pathWithChIP}/H3K27Ac_FL_E10.5+DMSO+6h.bw
title = H3K27ac
nans_to_zeros = true
min_value = 0
max_value = ${ymax_values_k27[$i]}
color = cyan
height = 3

[x-axis]
fontsize = 12
" >> ${ini_file}
done

pyGenomeTracks --tracks Axin2_Supplement.ini --region chr11:108350000-109350000 -o Axin2_Supplement.pdf

pyGenomeTracks --tracks Grem1_Supplement.ini --region chr2:113400000-113830000 -o Grem1_Supplement.pdf

pyGenomeTracks --tracks Bmp4_Supplement.ini --region chr14:45700000-46900000 -o Bmp4_Supplement.pdf

pyGenomeTracks --tracks Fgf8_Supplement.ini --region chr19:45425000-46000000 -o Fgf8_Supplement.pdf

pyGenomeTracks --tracks Shh_Supplement.ini --region chr5:28280000-29500000 -o Shh_Supplement.pdf  
