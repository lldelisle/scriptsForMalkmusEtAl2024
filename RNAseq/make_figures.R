# Load dependencies
if (!"devtools" %in% installed.packages()) {
    install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("pheatmap")
safelyLoadAPackageInCRANorBioconductor("RColorBrewer")
safelyLoadAPackageInCRANorBioconductor("reshape2")
safelyLoadAPackageInCRANorBioconductor("biomaRt")
safelyLoadAPackageInCRANorBioconductor("ggrepel")
safelyLoadAPackageInCRANorBioconductor("goseq")
safelyLoadAPackageInCRANorBioconductor("TxDb.Mmusculus.UCSC.mm10.ensGene")
safelyLoadAPackageInCRANorBioconductor("ggpubr")

# Variables:
gitHubDirectory <- "/home/ldelisle/Documents/mygit/scriptsForMalkmusEtAl2024/"
table.with.norm.exp.fn <- file.path(gitHubDirectory, "RNAseq", "outputs", "AllCufflinks_Simplified_timecourse.txt.gz")
samples.plan.fn <- file.path(gitHubDirectory, "RNAseq", "samplesPlan_time_course.txt")
list.colors <- list(
    "treatment" = c("DMSO" = "darkgreen", "C59" = "orange"),
    "timeAfterInjection" = colorRampPalette(c("cyan", "magenta"))(5),
    "cluster" = c("A" = "darkblue", "B" = "red", "C" = "lightblue")
)
names(list.colors$timeAfterInjection) <- paste0(c(1, 6, 12, 18, 24), "h")
summary.de.fn <- file.path(gitHubDirectory, "RNAseq", "outputs", "summary_significant.txt.gz")
l2fc.threshold <- 1
gene.panel.fn <- file.path(gitHubDirectory, "RNAseq", "gene_panel_list.txt")
n.var.genes <- 2000
output.directory <- file.path(gitHubDirectory, "RNAseq", "outputs")


dir.create(output.directory, showWarnings = FALSE, recursive = TRUE)

# Prepare inputs
## Samples plan
samples.plan.df <- read.delim(samples.plan.fn)
rownames(samples.plan.df) <- samples.plan.df$sample
# Reorder timeAfterInjection
samples.plan.df$timeAfterInjection <- factor(samples.plan.df$timeAfterInjection,
    levels = unique(samples.plan.df$timeAfterInjection)[order(as.numeric(gsub(
        "h", "",
        unique(samples.plan.df$timeAfterInjection)
    )))]
)
# Reorder treatment
samples.plan.df$treatment <- factor(samples.plan.df$treatment, levels = unique(samples.plan.df$treatment))

## Expression
expression.df <- read.delim(table.with.norm.exp.fn, check.names = FALSE)
# Exclude chrX, chrY, chrM
expression.df <- subset(expression.df, !grepl("chr[XYM]:", locus))
expression.mat <- expression.df[, paste0("FPKM_", rownames(samples.plan.df))]
rownames(expression.mat) <- expression.df$gene_id
colnames(expression.mat) <- gsub("^FPKM_", "", colnames(expression.mat))
expression.mat.l2 <- log2(1 + expression.mat)

## Annotations for genes
# Summary.df is with l2fc threshold of 1
summary.df <- read.delim(summary.de.fn)

# Make PCA and clustering for time series samples
current.sp <- subset(samples.plan.df, study == "Time-series")
current.l2.exp.mat <- expression.mat.l2[, rownames(current.sp)]
if (n.var.genes == "all") {
    current.l2.exp.mat.subset <- current.l2.exp.mat
} else {
    current.l2.exp.mat.subset <- current.l2.exp.mat[order(apply(current.l2.exp.mat, 1, var),
        decreasing = T
    )[1:min(nrow(current.l2.exp.mat), n.var.genes)], ]
}
sample.pca <- prcomp(t(current.l2.exp.mat.subset),
    center = TRUE,
    scale. = FALSE
)


sample.pca.meta <- data.frame(current.sp, sample.pca$x[rownames(current.sp), ])
var <- round((sample.pca$sdev)^2 / sum(sample.pca$sdev^2) * 100)
for (i in 1) {
    ggplot(sample.pca.meta, aes_string(x = paste0("PC", i), y = paste0("PC", i + 1))) +
        geom_point(aes_string(color = "treatment", shape = "timeAfterInjection"), size = 3) +
        xlab(paste0("PC", i, ":", var[i], "% variance")) +
        ylab(paste0("PC", i + 1, ":", var[i + 1], "% variance")) +
        scale_color_manual(values = list.colors[["treatment"]]) +
        theme_classic()
    ggsave(file.path(output.directory, paste0("PC", i, "_", i + 1, "_", n.var.genes, "vargenes.pdf")),
        width = 4, height = 3
    )
}
correlation.matrix <- cor(current.l2.exp.mat.subset, method = "spearman")
rownames(correlation.matrix) <- simplifiedNamesByStartEnd(rownames(correlation.matrix))
colnames(correlation.matrix) <- simplifiedNamesByStartEnd(colnames(correlation.matrix))
sample.dist <- as.dist(1 - correlation.matrix)
# Reorder clustering
clu <- hclust(sample.dist, method = "ward.D2")
dd <- as.dendrogram(clu)
clu2 <- reorder(dd, 1:nrow(correlation.matrix), agglo.FUN = min)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
annot <- current.sp[, c("treatment", "timeAfterInjection")]
rownames(annot) <- simplifiedNamesByStartEnd(rownames(annot))

pheatmap(
    correlation.matrix,
    cluster_rows = as.hclust(clu2),
    cluster_cols = as.hclust(clu2),
    cellwidth = 10,
    cellheight = 10,
    annotation_row = annot,
    annotation_col = annot,
    annotation_colors = list.colors,
    main = "spearmanCor - ward clustering",
    clustering_method = "ward.D2",
    col = rev(colors),
    filename = file.path(output.directory, paste0("spearman_cor_on", n.var.genes, ".pdf"))
)

# Count the number of DE genes in time-series
pairs <- unique(samples.plan.df$pair_group[samples.plan.df$study == "Time-series"])
nb.genes <- NULL
for (pair in pairs) {
    nb.up <- sum(!is.na(summary.df[, paste0(pair, "_padj")]) & summary.df[, paste0(pair, "_padj")] < 0.05 & summary.df[, paste0(pair, "_l2fc")] > l2fc.threshold)
    nb.down <- sum(!is.na(summary.df[, paste0(pair, "_padj")]) & summary.df[, paste0(pair, "_padj")] < 0.05 & summary.df[, paste0(pair, "_l2fc")] < -l2fc.threshold)
    nb.genes <- rbind(nb.genes, c(pair, "C59", nb.up), c(pair, "DMSO", nb.down))
}
nb.genes <- as.data.frame(nb.genes)
colnames(nb.genes) <- c("pair", "regulation", "nb")
nb.genes$nb <- as.numeric(nb.genes$nb)
nb.genes$pair <- factor(nb.genes$pair, levels = unique(nb.genes$pair))
nb.genes$timeAfterInjection <- gsub("DMSOvsC59_", "", nb.genes$pair)
nb.genes$timeAfterInjection <- factor(nb.genes$timeAfterInjection, levels = levels(samples.plan.df$timeAfterInjection))
print(nb.genes)
#             pair regulation  nb timeAfterInjection
# 1   DMSOvsC59_1h        C59  16                 1h
# 2   DMSOvsC59_1h       DMSO  13                 1h
# 3   DMSOvsC59_6h        C59  37                 6h
# 4   DMSOvsC59_6h       DMSO  73                 6h
# 5  DMSOvsC59_12h        C59 428                12h
# 6  DMSOvsC59_12h       DMSO 203                12h
# 7  DMSOvsC59_18h        C59  78                18h
# 8  DMSOvsC59_18h       DMSO 107                18h
# 9  DMSOvsC59_24h        C59  54                24h
# 10 DMSOvsC59_24h       DMSO  51                24h
ggplot(nb.genes, aes(x = timeAfterInjection, y = nb)) +
    geom_bar(aes(fill = regulation), stat = "identity") +
    ylab("Number of DE genes") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    xlab("Time after injection") +
    scale_fill_manual("increased in", values = list.colors$treatment)
ggsave(file.path(output.directory, "nb_genes.pdf"), width = 4, height = 4)

# Study evolution of l2fc across time of treatment:
# I use genes which are significant in at least one time
my.genes.id <- unique(unlist(lapply(pairs, function(pair) {
    summary.df$gene_id[!is.na(summary.df[, paste0(pair, "_padj")]) & summary.df[, paste0(pair, "_padj")] < 0.05 & abs(summary.df[, paste0(pair, "_l2fc")]) > l2fc.threshold]
})))
my.genes <- summary.df$gene_short_name[match(my.genes.id, summary.df$gene_id)]

# Genes downregulated at 6 hours
pair <- "DMSOvsC59_6h"
cat(sort(summary.df$gene_short_name[!is.na(summary.df[, paste0(pair, "_padj")]) & summary.df[, paste0(pair, "_padj")] < 0.05 & summary.df[, paste0(pair, "_l2fc")] < -l2fc.threshold]), sep = ", ")
cat("\n")

# Extract log2fc
mat.all <- summary.df[, grep("^DMSOvsC59_.*h_l2fc", colnames(summary.df))]
mat <- mat.all[match(my.genes, summary.df$gene_short_name), ]
rownames(mat) <- my.genes
# All values with NA means it was not expressed in any of the conditions:
mat[is.na(mat)] <- 0

# Cluster genes by log2fc values:
meth <- "ward.D2"
nClusters <- 3
sub.df.scale <- t(apply(mat, 1, scale))
colnames(sub.df.scale) <- colnames(mat)
d <- as.dist(1 - cor(t(sub.df.scale)))
hc <- hclust(d, method = meth)
annot.genes.cluster <- data.frame(ori.cluster = LETTERS[1:nClusters][cutree(hc, nClusters)])
rownames(annot.genes.cluster) <- hc$labels
# We switch 'A' and 'B'
new.label <- c("B" = "A", "A" = "B", "C" = "C")
annot.genes.cluster$cluster <- new.label[annot.genes.cluster$ori.cluster]
annot.genes.cluster$ori.cluster <- NULL

# Plot a heatmap
pheatmap(sub.df.scale,
    cluster_cols = F,
    color = colorRampPalette(c("navy", "white", "red"))(50),
    breaks = c(min(-2.01, min(apply(mat, 1, scale))), seq(-2, 2, length.out = 48), max(2.01, max(apply(mat, 1, scale)))),
    annotation_row = annot.genes.cluster,
    annotation_colors = list.colors,
    show_rownames = F,
    clustering_distance_rows = "correlation",
    clustering_method = meth,
    main = "Log2FC centered scaled clustered by correlation\nAll genes significant at one time point",
    file = file.path(output.directory, paste0("Timeseries_l2fc_cs_clustering_cor_", meth, "_", nClusters, "clust.pdf")),
    annotation_names_row = FALSE
)

# Show the log2fc cluster per cluster
df.gg <- melt(sub.df.scale)
df.gg.non.scale <- melt(as.matrix(mat))
colnames(df.gg.non.scale)[3] <- "l2fc"
df.gg <- merge(df.gg, df.gg.non.scale)
df.gg$time <- as.numeric(gsub("^DMSOvsC59_", "", gsub("h_l2fc", "", df.gg$Var2)))
df.gg$cluster <- annot.genes.cluster[as.character(df.gg$Var1), "cluster"]
ggplots <- list()
for (my.cluster in unique(df.gg$cluster)) {
    ggplots[[my.cluster]] <-
        ggplot(subset(df.gg, cluster == my.cluster), aes(x = time, y = l2fc)) +
        geom_line(aes(group = Var1), alpha = 0.01) +
        stat_summary(aes(color = cluster), fun = mean, geom = "line", show.legend = FALSE) +
        stat_summary(aes(color = cluster), fun = mean, geom = "point", show.legend = FALSE) +
        theme_minimal() +
        coord_cartesian(ylim = c(-2, 2)) +
        scale_color_manual("", values = list.colors[["cluster"]]) +
        scale_x_continuous(name = "time [h]", breaks = c(1, 6, 12, 18, 24), minor_breaks = NULL) +
        ylab("log2(C59/DMSO)")
}

# GO on each cluster
# We use the same version as gtf (102):
mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "https://nov2020.archive.ensembl.org")
mart <- useDataset("mmusculus_gene_ensembl", mart)
genes <- getBM(
    attributes = c("ensembl_gene_id", "gene_biotype", "external_gene_name", "go_id", "namespace_1003"),
    mart = mart
)

# GO per cluster and store to file
if (!file.exists(file.path(output.directory, paste0("k", nClusters, "_cluster_", letters[nClusters], "_BP.txt")))) {
    pc.genes <- subset(genes, gene_biotype == "protein_coding")
    genes.length <- getlength(summary.df$gene_id, "mm10", "ensGene")
    gene2cat.bp <- with(
        subset(pc.genes, namespace_1003 == "biological_process"),
        split(go_id, ensembl_gene_id)
    )
    for (my.cluster in unique(annot.genes.cluster[, "cluster"])) {
        signif.genes <- rownames(subset(annot.genes.cluster, annot.genes.cluster[, "cluster"] == my.cluster))
        gene.vector <- as.integer(summary.df$gene_short_name %in% signif.genes)
        names(gene.vector) <- summary.df$gene_id
        pwf <- nullp(gene.vector, bias.data = genes.length)
        GO.BP <- goseq(pwf, gene2cat = gene2cat.bp)
        write.table(GO.BP, file.path(output.directory, paste0("k", nClusters, "_cluster_", my.cluster, "_BP.txt")), sep = "\t", quote = F, row.names = F)
        dev.off()
    }
}

# Plot the GO result
for (my.cluster in unique(annot.genes.cluster$cluster)) {
    GO.BP <- read.delim(file.path(output.directory, paste0("k", nClusters, "_cluster_", my.cluster, "_BP.txt")))
    GO.BP$term <- factor(GO.BP$term, levels = rev(unique(GO.BP$term)))
    n <- 9
    ggplots[[paste0(my.cluster, "_GO")]] <- ggplot(GO.BP[1:n, ], aes(x = -log10(over_represented_pvalue), y = term)) +
        geom_col(fill = "#6395ecff", color = "#5d69e1ff") +
        xlab("-log10(p-value)") +
        ylab("") +
        scale_y_discrete(position = "right") +
        scale_x_continuous(expand = c(0, 0)) +
        theme_minimal()
}
g <- ggarrange(
    plotlist = ggplots[c("A", "A_GO", "B", "B_GO", "C", "C_GO")],
    ncol = 2, nrow = 3, widths = c(1, 2)
)
ggsave(file.path(output.directory, "profile_GO.pdf"), width = 8, height = 8)

# Display selected list:
gene.panel.df <- read.delim(gene.panel.fn)
# setdiff(gene.panel.df$gene, summary.df$gene_short_name)
input.mat <- mat.all[match(unique(gene.panel.df$gene), summary.df$gene_short_name), ]
rownames(input.mat) <- unique(gene.panel.df$gene)
input.mat[is.na(input.mat)] <- 0
df.gg.non.scale.selected <- melt(as.matrix(input.mat))
colnames(df.gg.non.scale.selected) <- c("gene", "Var2", "l2fc")
df.gg.non.scale.selected$time <- as.numeric(gsub("^DMSOvsC59_", "", gsub("h_l2fc", "", df.gg.non.scale.selected$Var2)))
df.gg.non.scale.selected$gene <- as.character(df.gg.non.scale.selected$gene)

for (my.panel in unique(gene.panel.df$panel)) {
    input.df <- merge(
        df.gg.non.scale.selected,
        subset(
            gene.panel.df,
            panel == my.panel
        )
    )
    input.df$label <- ifelse(
        input.df$highlight == "yes",
        input.df$gene,
        "other"
    )
    input.df$label <- factor(input.df$label, levels = c(setdiff(sort(unique(input.df$label)), "other"), "other"))
    colors <- c(scales::hue_pal()(length(levels(input.df$label)) - 1), "grey70")
    names(colors) <- levels(input.df$label)
    avg.18 <- mean(input.df$l2fc[input.df$time == 18])
    if (avg.18 > -1) {
        y <- 0
    } else {
        y <- 1
    }
    if ("other" %in% input.df$label) {
        linewidth.legend <- c(rep(0.5, nlevels(input.df$label) - 1), 0.25)
    } else {
        linewidth.legend <- 0.5
    }
    g <- ggplot(
        input.df[rev(order(input.df$label)), ],
        aes(x = time, y = l2fc, color = label)
    ) +
        geom_line(aes(group = gene, linewidth = label == "other")) +
        theme_minimal() +
        coord_cartesian(ylim = c(-4, 2)) +
        # scale_color_manual("", values = list.colors[["cluster"]]) +
        scale_x_continuous(name = "time [h]", breaks = c(1, 6, 12, 18, 24), minor_breaks = NULL) +
        ylab("log2(C59/DMSO)") +
        ggtitle(my.panel) +
        guides(
            linewidth = "none",
            colour = guide_legend(
                title = "",
                override.aes = list(linewidth = linewidth.legend)
            )
        ) +
        scale_color_manual(values = colors) +
        scale_linewidth_manual(values = c("TRUE" = 0.25, "FALSE" = 0.5))
    if (my.panel == "fgf and HH main") {
        g1 <- g +
                coord_cartesian(ylim = c(-3.5, 2)) +
                theme(
                    axis.text.x = element_blank(),
                    axis.title.x = element_blank()
                )
        g2 <- g +
            coord_cartesian(ylim = c(-6.5, -5.5)) +
            scale_y_continuous(breaks = -6) +
            ylab("") +
            ggtitle("") +
            theme(legend.position = "none")
        # build the plots
        p1.common.x <- ggplot_gtable(ggplot_build(g1))
        p2.common.x <- ggplot_gtable(ggplot_build(g2))

        # copy the plot width from p1 to p2
        p2.common.x$widths <- p1.common.x$widths

        g <- ggarrange(
            p1.common.x, p2.common.x,
            nrow = 2, heights = c(1, 0.3)
        )
    }
    ggsave(file.path(output.directory, paste0("line_plot_", gsub(" ", "_", my.panel), ".pdf")),
    g,
    width = 4.5, height = 4
)
}

# Export the annotations
interesting.categories <- c(
    "Cartilage_GO:0051216", "HH_GO:0007224", "WNT_GO:0016055", "BMP_GO:0030509", "TGFb_GO:0007179", "FGFR_GO:0008543"
)

annot.df <- summary.df[summary.df$gene_short_name %in% rownames(annot.genes.cluster), 1:2]
annot.df$cluster <- annot.genes.cluster[annot.df$gene_short_name, "cluster"]
for (ic in interesting.categories) {
    go.term <- gsub("^.*_", "", ic)
    annot.df[, ic] <- ifelse(
        annot.df$gene_id %in% genes$ensembl_gene_id[genes$go_id == go.term],
        "yes", "no"
    )
}
annot.df$panel <- gene.panel.df$panel[match(annot.df$gene_short_name, gene.panel.df$gene)]
annot.df$highlight <- gene.panel.df$highlight[match(annot.df$gene_short_name, gene.panel.df$gene)]
write.table(annot.df, "RNAseq/gene_annotations.txt", sep = "\t", quote = FALSE, row.names = FALSE)
