if (!"devtools" %in% installed.packages()) {
    install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("GenomicRanges")
safelyLoadAPackageInCRANorBioconductor("rtracklayer")
safelyLoadAPackageInCRANorBioconductor("combinat")
safelyLoadAPackageInCRANorBioconductor("BiocGenerics")
devtools::install_github("lldelisle/usefulLDfunctionsGR")
library(usefulLDfunctionsGR)

# Variables
gitHubDirectory <- "/home/ldelisle/Documents/mygit/scriptsForMalkmusEtAl2024/"
url.gtf <- "https://zenodo.org/records/7510406/files/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz?download=1"
gtf.path <- "mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz"
url.blacklist <- "https://github.com/Boyle-Lab/Blacklist/raw/f4a45ab5b52d427a8689a1a143e8504e1cba5962/lists/mm10-blacklist.v2.bed.gz"
black.list.path <- file.path(gitHubDirectory, "ChIPseq/outputs/mm10-blacklist.v2.bed.gz")
vista.path <- file.path(gitHubDirectory, "ChIPseq/outputs/limb_vista_mouse_mm10.bed")
output.directory <- file.path(gitHubDirectory, "ChIPseq/outputs/")
peaks.directory <- commandArgs(TRUE)[1]
dist.to.TSS <- 1000
dist.to.Ac <- 3000
my.sample <- "FL_E10.5"

# Get inputs
if (!file.exists(gtf.path)) {
    download.file(url.gtf, destfile = gtf.path, mode = "wb")
}
tss <- getTSSinUCSCFormatFromEnsemblGTF(gtf.path)
if (!file.exists(black.list.path)) {
    download.file(url.blacklist, destfile = black.list.path, mode = "wb")
}
black.list <- import.bed(black.list.path)
vista.list <- import.bed(vista.path)
bcat.chip.peak.fn <- file.path(gitHubDirectory, "ChIPseq", "outputs", paste0("b-Catenin_", my.sample, "+DMSO+6h_consensus.narrowPeak.gz"))
# From https://charlesjb.github.io/How_to_import_narrowPeak/
extraCols_narrowPeak <- c(
    signalValue = "numeric",
    pValue = "numeric",
    qValue = "numeric",
    relSummit = "integer"
)
gr_narrowPeak <- import(
    bcat.chip.peak.fn,
    format = "BED",
    extraCols = extraCols_narrowPeak
)
bcat.chip.summits <- GenomicRanges::resize(
    shift(
        gr_narrowPeak,
        shift = gr_narrowPeak$relSummit
    ),
    width = 1, fix = "start"
)
gr_narrowPeak_Ac1 <- import(
    file.path(peaks.directory, paste0("H3K27Ac_", my.sample, "+DMSO+6h_rep1.narrowPeak.gz")),
    format = "BED",
    extraCols = extraCols_narrowPeak
)
gr_narrowPeak_Ac2 <- import(
    file.path(peaks.directory, paste0("H3K27Ac_", my.sample, "+DMSO+6h_rep2.narrowPeak.gz")),
    format = "BED",
    extraCols = extraCols_narrowPeak
)
# Check intersection with Vista
mcols(vista.list)[, paste0(my.sample, "_intersect")] <-
    overlapsAny(
        vista.list, bcat.chip.summits
    )
cat("Remove", sum(overlapsAny(bcat.chip.summits, black.list)), "summits because they overlap black list.\n")
bcat.chip.summits <- subsetByOverlaps(
    bcat.chip.summits, black.list,
    invert = TRUE
)
mcols(vista.list)[, paste0(my.sample, "_intersect_noBlacklist")] <-
    overlapsAny(
        vista.list, bcat.chip.summits
    )
bcat.chip.summits$TSS_intersect <-
    overlapsAny(
        bcat.chip.summits, tss,
        maxgap = dist.to.TSS
    )
bcat.chip.summits$VistaLimb_intersect <-
    overlapsAny(
        bcat.chip.summits, vista.list
    )
bcat.chip.summits$H3K27Ac_intersect <-
    overlapsAny(
        bcat.chip.summits, gr_narrowPeak_Ac1,
        maxgap = dist.to.Ac
    ) | overlapsAny(
        bcat.chip.summits, gr_narrowPeak_Ac2,
        maxgap = dist.to.Ac
    )
bcat.chip.summits$category <-
    apply(
        mcols(bcat.chip.summits)[, c("TSS_intersect", "VistaLimb_intersect", "H3K27Ac_intersect")],
        1,
        function(v) {
            gsub("_intersect", "", paste(names(v)[v], collapse = "&"))
        }
    )
print(table(bcat.chip.summits$category))
export.bed(
    subset(bcat.chip.summits, TSS_intersect),
    con = file.path(output.directory, paste0("bcat_with_TSS_", my.sample, ".bed"))
)
export.bed(
    subset(bcat.chip.summits, !TSS_intersect),
    con = file.path(output.directory, paste0("bcat_noTSS_", my.sample, ".bed"))
)
write.table(
    as.data.frame(bcat.chip.summits),
    file.path(output.directory, paste0("bcat_", my.sample, "_with_annotations.txt")),
    quote = FALSE, sep = "\t", row.names = FALSE
)
# Remove 253 summits because they overlap black list.
#                                 H3K27Ac                   TSS 
#              8227                  5188                   232 
#       TSS&H3K27Ac TSS&VistaLimb&H3K27Ac             VistaLimb 
#              2147                    15                    26 
# VistaLimb&H3K27Ac 
#               150 
write.table(
    as.data.frame(vista.list),
    file.path(output.directory, paste0("limb_vista_with_annotations.txt")),
    quote = FALSE, sep = "\t", row.names = FALSE
)
