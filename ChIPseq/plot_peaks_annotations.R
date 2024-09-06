library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("ggplot2")
my.colors <- c(
    "no-H3K27Ac" = "grey30",
    "H3K27Ac" = "cyan",
    "no-bCat" = "grey30",
    "bCat" = "red"
)
output.directory <- "ChIPseq/outputs"
my.sample <- "FL_E10.5"
df.annot <- read.delim(file.path(output.directory, paste0("bcat_", my.sample, "_with_annotations.txt")))
df.annot$TSS_intersect_pretty <- ifelse(
    df.annot$TSS_intersect,
    "< 1kb",
    "> 1kb"
)
df.annot$H3K27Ac_intersect_pretty <- ifelse(
    df.annot$H3K27Ac_intersect,
    "H3K27Ac",
    "no-H3K27Ac"
)
ggplot(df.annot, aes(x = TSS_intersect_pretty, fill = H3K27Ac_intersect_pretty)) +
    geom_bar(position = position_dodge()) +
    geom_text(
        aes(label = after_stat(count), y = after_stat(count)),
        stat = "count", position = position_dodge(0.9),
        vjust = -.7
    ) +
    ylab(expression("Number of " * beta * "-Catenin summits")) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    xlab("from TSS") +
    scale_fill_manual("Overlap", values = my.colors) +
    theme_classic()
ggsave(file.path(output.directory, paste0("bcat_", my.sample, "_annotations.pdf")), width = 4, height = 4)
print(with(df.annot, table(TSS_intersect_pretty, H3K27Ac_intersect_pretty)))
#                     H3K27Ac_intersect_pretty
# TSS_intersect_pretty H3K27Ac no-H3K27Ac
#                < 1kb    2162        232
#                > 1kb    5338       8253

# Vista proportion
df.annot <- read.delim(file.path(
    output.directory,
    "limb_vista_with_annotations.txt"
))
df.annot[, paste0(my.sample, "_intersect_pretty")] <- ifelse(
    df.annot[, paste0(my.sample, "_intersect")],
    "bCat",
    "no-bCat"
)
ggplot(
    df.annot,
    aes(
        x = .data[[paste0(my.sample, "_intersect_pretty")]],
        fill = .data[[paste0(my.sample, "_intersect_pretty")]]
    )
) +
    geom_bar(show.legend = FALSE) +
    geom_text(
        aes(label = after_stat(count), y = after_stat(count)),
        stat = "count",
        vjust = -.7
    ) +
    geom_text(
        aes(label = paste0("(", scales::percent(after_stat(count) / nrow(df.annot)), ")"), y = after_stat(count)),
        stat = "count",
        vjust = 1.4
    ) +
    ylab("Number of Limb enhancer (VISTA)") +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    xlab("") +
    scale_fill_manual("", values = my.colors) +
    theme_classic() +
    ggtitle(my.sample)
ggsave(file.path(output.directory, paste0("limb_vista_", my.sample, "_annotations.pdf")), width = 4, height = 4)
