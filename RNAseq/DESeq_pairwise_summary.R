options(stringsAsFactors = F)
rm(list = ls())

# Get arguments
gitHubDirectory <- commandArgs(TRUE)[1]
path <- commandArgs(TRUE)[2]
# Load dependencies
if (!"devtools" %in% installed.packages()) {
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("DESeq2")
safelyLoadAPackageInCRANorBioconductor("rtracklayer")
# Wrapper for DESeq2
simpleDeseqAna <- function(count.table, factorForAna, pathOutput,
                           samplesPlan, LRT=F, pvalT=0.05,
                           lfcT=1.5, writeRLOG=F, ...){
  # Checking the conditions
  if (!(factorForAna %in% colnames(samplesPlan))) {
    stop("This is not part of the column names.")
  }
  if (length(levels(samplesPlan[, factorForAna])) == 1) {
    stop("The factor you chose have only 1 value. The analysis is not possible.")
  }
  if (length(levels(samplesPlan[, factorForAna])) > 2 & !LRT) {
    print("The factor you chose have more than 2 values. LRT will be applied.")
    LRT <- T
  }
  # Lauching DESeq2
  cmd <- paste0("dds <- DESeqDataSetFromMatrix(countData = count.table[, match(rownames(samplesPlan), colnames(count.table))],
             colData = samplesPlan,
             design = ~ ", factorForAna, ")")
  eval(parse(text = cmd))
  print("Genes that are never expressed are removed")
  dds <- dds[rowSums(counts(dds)) > 1, ]
  if (LRT) {
    dds <- DESeq(dds, minReplicatesForReplace = Inf, test = "LRT", reduced = ~ 1)
  } else{
    dds <- DESeq(dds, minReplicatesForReplace = Inf)
  }
  res <- results(dds, ...)
  resOrdered <- res[order(res$padj), ]
  # Subsetting the annotation file
  ann <- count.table[, colnames(count.table) %in% c("gene_id", "gene_short_name", "locus")]
  annot.df <- data.frame(ann[match(rownames(resOrdered), ann$gene_id), ])
  resToExport <- data.frame(annot.df, counts(dds, normalized = TRUE)[rownames(resOrdered), ], resOrdered)
  write.table(resToExport, file = paste0(pathOutput, "DESeq2Results.txt"), sep = '\t', row.names = F, quote = F)
  dfDiffExp <- subset(resToExport, resToExport$padj < pvalT & abs(resToExport$log2FoldChange) > lfcT)
  write.table(dfDiffExp, file = paste0(pathOutput, "DESeq2significant.txt"), sep = '\t', row.names = F, quote = F)
  rld <- rlog(dds)
  rlogdata <- assay(rld)
  if (writeRLOG) {
    resToExport2 <- data.frame(annot.df, rlogdata[rownames(resOrdered), ])
    write.table(resToExport2, file = paste0(pathOutput, "rlog.txt"), sep = '\t', row.names = F, quote = F)
  }
  return(invisible(dfDiffExp))
}

# Fixed variables:
tableWithCounts <- file.path(path, "mergedTables/AllHTSeqCounts_subset.txt")
gtf.file <- list.files(path, pattern = "merge.*UCSC.gtf", full.names = T)
samplesPlan <- file.path(gitHubDirectory, "RNAseq", "samplesPlan_time_course.txt")
factorToStudy <- "treatmentOrGenotype"
pathForDESeq2 <- file.path(path, "DESeq2_pairwise/")
log2FC.threshold <- log2(2)

# Prepare inputs
samplesPlanDF <- read.delim(samplesPlan)
rownames(samplesPlanDF) <- samplesPlanDF$sample
factorizedSP <- samplesPlanDF
for (cn in colnames(factorizedSP)) {
  uniqVal <- unique(factorizedSP[, cn])
  factorizedSP[, cn] <- factor(factorizedSP[, cn], levels = uniqVal)
}
samplesToPlot <- samplesPlanDF$sample
count.table <- read.delim(tableWithCounts, check.names = FALSE)
colnames(count.table)[1] <- "gene_id"

if (!dir.exists(pathForDESeq2)) {
  dir.create(pathForDESeq2, recursive = T)
}
# Prepare a big table with the results of all DESeq2
gtf <- readGFF(gtf.file)
big.annot <- unique(gtf[, c("gene_id", "gene_name", "seqid", "gene_biotype")])
rm(gtf)
colnames(big.annot) <- c("gene_id", "gene_short_name", "chr", "gene_biotype")
count.table <- merge(count.table, big.annot, all.x = T)
# I keep only protein_coding genes
count.table <- subset(count.table, gene_biotype == "protein_coding")
rownames(count.table) <- count.table$gene_id
big.annot <- subset(big.annot, gene_biotype == "protein_coding")

for (study in unique(samplesPlanDF$study)) {
  print(study)
  big.annot[, study] <- ""
  for (pair in unique(samplesPlanDF$pair_group[samplesPlanDF$study == study])) {
    print(pair)
    # Select the samples
    new.samples.plan <- subset(samplesPlanDF, pair_group == pair)
    if ("DMSO" %in% new.samples.plan[, factorToStudy]) {
      new.samples.plan[, factorToStudy] <- factor(new.samples.plan[, factorToStudy],
                                                  levels = c("DMSO", setdiff(unique(new.samples.plan[, factorToStudy]), "DMSO")))
    }
    # Run or read DESeq2 results with Wald test threshold of FC specified above
    if ( !file.exists(paste0(pathForDESeq2, "/", study, "_", pair, "_DESeq2significant.txt"))) {
      signif <- simpleDeseqAna(count.table, factorToStudy, paste0(pathForDESeq2,"/", study, "_", pair, "_"),
                               new.samples.plan,
                               LRT = F, lfcT = log2FC.threshold, writeRLOG = F) #,
                               # theta = c(0.15, 0.99))
    } else {
      signif <- read.delim(paste0(pathForDESeq2, "/", study, "_", pair, "_DESeq2significant.txt"))
    }
    full.res <- read.delim(paste0(pathForDESeq2, "/", study, "_", pair, "_DESeq2Results.txt"))
    # Add results to the dataframe
    big.annot[, paste0(pair, "_l2fc")] <- full.res$log2FoldChange[match(big.annot$gene_id, full.res$gene_id)]
    big.annot[, paste0(pair, "_padj")] <- full.res$padj[match(big.annot$gene_id, full.res$gene_id)]
    big.annot[big.annot$gene_id %in% signif$gene_id, study] <- paste0(big.annot[big.annot$gene_id %in% signif$gene_id, study], "_", pair)
  }
  big.annot[, study] <- gsub("^_", "", big.annot[, study])
}
write.table(big.annot, file.path(pathForDESeq2, "summary_significant.txt"), sep = "\t", quote = F, row.names = F)
