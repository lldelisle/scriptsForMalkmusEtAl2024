# Get vista enhancers in limb:
mkdir -p ChIPseq/outputs
wget "https://enhancer.lbl.gov/cgi-bin/imagedb3.pl?page=1;search.expression=26;search.status=Positives;search.result=yes;show=1;page_size=100;form=ext_search;search.gene=;action=search;search.org=Mouse;search.sequence=1" -O "ChIPseq/outputs/limb_vista_mouse.fa"
# Remove the <pre>
sed -i -e "s/<pre>//g" -e "s#</pre>##" ChIPseq/outputs/limb_vista_mouse.fa
# Convert fasta to bed
grep ">" ChIPseq/outputs/limb_vista_mouse.fa | awk -F "|" '{split($2, a, ":|-"); gsub(" ", "", $3); print a[1]"\t"a[2]"\t"a[3]"\t"$3}' > ChIPseq/outputs/limb_vista_mouse_mm9.bed
# Get liftOver for linux
wget "https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver" -O "/tmp/liftOver" -nc
chmod +x /tmp/liftOver
# Get the chainFile
wget "https://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz" -O /tmp/mm9ToMm10.over.chain.gz -nc
# liftOver
/tmp/liftOver ChIPseq/outputs/limb_vista_mouse_mm9.bed /tmp/mm9ToMm10.over.chain.gz ChIPseq/outputs/limb_vista_mouse_mm10.bed ChIPseq/outputs/limb_vista_mouse_mm9_unmapped.bed

# Intersect consensus peaks for b-cat with TSS, H3K27Ac and VISTA enhancers
Rscript ChIPseq/compare_peaks_annotations.R ~/mnt/scratch/Sheth/ChIPseq/toGEO/
