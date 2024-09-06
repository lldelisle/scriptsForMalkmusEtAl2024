# Inputs
# bcat_${sample}_noblacklist = cat bcat_noTSS_${sample}.bed bcat_with_TSS_${sample}.bed
# bcat_noTSS_${sample} = bcat_noTSS_${sample}.bed
# bcat_${sample}_noTSS_H3K27Ac = filter bcat_${sample}_with_annotations.txt for TRUE in column 14 and FALSE in column 12 and cut 1,2,3,4.

# Homer version 4.11 in galaxy
# Command line is
findMotifsGenome.pl 'bcat_FL_E10_5_noblacklist' 'mm10_UCSC.fa' 'bcat_FL_E10_5_noblacklist_motif' -preparsedDir '/data/galaxy/galaxy/var/tool-data/mm10_UCSC/homer_preparse/mm10_UCSC_200' -size '200' -len '8,10,12' -S 25 -mis 2    -mset auto    -noknown   -gc  -local 0 -redundant 2.0 -maxN 0.7   -nlen '3' -nmax '160'  -e '0.0'  -minlp '-10.0'