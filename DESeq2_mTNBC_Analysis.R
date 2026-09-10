

# Load packages
required_packages <- c(
  "DESeq2",
  "ComplexHeatmap",
  "circlize",
  "ggplot2",
  "ggrepel"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing packages: ",
    paste(missing_packages, collapse = ", "),
    "\nRun renv::restore() before running this script."
  )
}

library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggrepel)
library(grid)


# Project folders


dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# 1. Download and read the GEO count matrix
data_file <- file.path(
  "data",
  "GSE264108_readcounts.txt.gz"
)

if (!file.exists(data_file)) {
  message("Downloading GSE264108 count matrix from GEO...")
  download.file(
    url = paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/",
      "GSE264nnn/GSE264108/suppl/",
      "GSE264108_readcounts.txt.gz"
    ),
    destfile = data_file,
    mode = "wb"
  )
}
raw_data <- read.delim(
  gzfile(data_file),
  check.names = FALSE
)
rownames(raw_data) <- raw_data$ensembl_gene_id


# 2. Prepare the count matrix
# Columns 2-15 contain raw counts for the 14 samples
cts <- as.matrix(
  raw_data[, 2:15]
)

cat(
  "\nCount matrix:",
  nrow(cts),
  "genes x",
  ncol(cts),
  "samples\n"
)

print(colnames(cts))


# Basic checks before DESeq2
stopifnot(
  ncol(cts) == 14,
  !anyNA(cts),
  all(cts >= 0),
  all(cts == round(cts))
)


# The paper reports removing genes only when expression was zero across all samples
all_zero <- sum(
  rowSums(cts) == 0
)

cat(
  "\nGenes with all-zero expression:",
  all_zero,
  "\n"
)

cts <- cts[
  rowSums(cts) > 0,
  ,
  drop = FALSE
]


# 3. Create sample information

# Healthy donor samples contain "Healthy" in the deposited
# sample names. The remaining samples belong to the mTNBC group.
condition <- ifelse(
  grepl("Healthy", colnames(cts),ignore.case = TRUE),
  "Healthy",
  "mTNBC"
)

colData <- data.frame(condition = condition,
  row.names = colnames(cts)
)


# Healthy donors are the reference group
colData$condition <- factor(colData$condition,levels = c("Healthy","mTNBC"))
cat("\nSamples per group:\n")
print(table(colData$condition))


# This dataset should contain seven samples in each group
stopifnot(sum(colData$condition == "Healthy") == 7,sum(colData$condition == "mTNBC") == 7)

# Count columns and metadata must be in the same order
stopifnot(all(colnames(cts) ==rownames(colData)))

# 4. Differential-expression analysis
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = colData,
  design = ~ condition)


# DESeq2 works directly with the raw counts
dds <- DESeq(dds)


# mTNBC is compared with Healthy.
# Positive log2FoldChange = higher in mTNBC
# Negative log2FoldChange = higher in Healthy
res <- results(dds,
  contrast = c(
    "condition",
    "mTNBC",
    "Healthy"
  ),
  alpha = 0.05
)

# 5. Prepare the DESeq2 results
res_df <- as.data.frame(res)
res_df$ensembl_gene_id <- rownames(res_df)

# Add gene symbols from the GEO annotation
res_df$gene_symbol <- raw_data$external_gene_id[
  match(
    res_df$ensembl_gene_id,
    raw_data$ensembl_gene_id
  )
]


# Use gene symbols for display. If no symbol is available, keep the Ensembl ID.
res_df$gene_label <- res_df$gene_symbol

missing_label <- is.na(res_df$gene_label) |
  res_df$gene_label == ""

res_df$gene_label[missing_label] <-
  res_df$ensembl_gene_id[missing_label]


# Formal DEG definition used in this analysis
sig_df <- res_df[
  !is.na(res_df$padj) &
    res_df$padj < 0.05,
  ,
  drop = FALSE
]


cat(
  "\nDifferential-expression results\n",
  "--------------------------------\n",
  "Significant genes:",
  nrow(sig_df),
  "\nUpregulated:",
  sum(sig_df$log2FoldChange > 0),
  "\nDownregulated:",
  sum(sig_df$log2FoldChange < 0),
  "\n"
)


# Save the DESeq2 results
write.csv(
  res_df,
  file.path(
    "results",
    "DESeq2_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  sig_df,
  file.path(
    "results",
    "DESeq2_significant_genes.csv"
  ),
  row.names = FALSE
)



# 6. Automatically identify top significant DEGs

sig_up <- sig_df[
  sig_df$log2FoldChange > 0,
  ,
  drop = FALSE
]

sig_down <- sig_df[
  sig_df$log2FoldChange < 0,
  ,
  drop = FALSE
]


# Largest positive fold changes
sig_up <- sig_up[
  order(
    sig_up$log2FoldChange,
    decreasing = TRUE
  ),
]

sig_up <- sig_up[
  !duplicated(sig_up$gene_label),
]

top10_up <- head(
  sig_up,
  10
)


# Most negative fold changes
sig_down <- sig_down[
  order(
    sig_down$log2FoldChange,
    decreasing = FALSE
  ),
]

sig_down <- sig_down[
  !duplicated(sig_down$gene_label),
]

top10_down <- head(
  sig_down,
  10
)


cat(
  "\nTop 10 significant upregulated genes:\n"
)

print(
  top10_up[
    ,
    c(
      "gene_label",
      "log2FoldChange",
      "pvalue",
      "padj"
    )
  ]
)


cat(
  "\nTop 10 significant downregulated genes:\n"
)

print(
  top10_down[
    ,
    c(
      "gene_label",
      "log2FoldChange",
      "pvalue",
      "padj"
    )
  ]
)


write.csv(
  top10_up,
  file.path(
    "results",
    "top10_upregulated_genes.csv"
  ),
  row.names = FALSE
)

write.csv(
  top10_down,
  file.path(
    "results",
    "top10_downregulated_genes.csv"
  ),
  row.names = FALSE
)


# Figure 5a - Heatmap

# The paper used DESeq2 for differential-expression testing, but Qlucore for visualization.
# GEO provides a separate normalized count matrix, so we use
# that matrix here for the heatmap rather than changing the
# raw counts used by DESeq2.

normalized_file <- file.path(
  "data",
  "GSE264108_readcounts_normalized_10mil.txt.gz"
)

if (!file.exists(normalized_file)) {
  
  message("Downloading normalized count matrix from GEO...")
  
  download.file(
    url = paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/",
      "GSE264nnn/GSE264108/suppl/",
      "GSE264108_readcounts_normalized_10mil.txt.gz"
    ),
    destfile = normalized_file,
    mode = "wb"
  )
}


normalized_data <- read.delim(
  gzfile(normalized_file),
  check.names = FALSE
)

rownames(normalized_data) <-
  normalized_data$ensembl_gene_id


# The same 14 sample columns
normalized_mat <- as.matrix(
  normalized_data[, 2:15]
)


# Keep the significant genes identified by DESeq2
heatmap_mat <- normalized_mat[
  rownames(sig_df),
  ,
  drop = FALSE
]


# Healthy donors first, followed by mTNBC
healthy_samples <- colnames(heatmap_mat)[
  grepl(
    "Healthy",
    colnames(heatmap_mat),
    ignore.case = TRUE
  )
]

mtnbc_samples <- setdiff(
  colnames(heatmap_mat),
  healthy_samples
)


# Cluster samples within each group.
# This keeps HDs together and mTNBC together while still
# allowing the samples inside each group to cluster naturally.

hd_tree <- hclust(
  dist(
    t(heatmap_mat[, healthy_samples, drop = FALSE])
  ),
  method = "average"
)

mtnbc_tree <- hclust(
  dist(
    t(heatmap_mat[, mtnbc_samples, drop = FALSE])
  ),
  method = "average"
)


healthy_samples <- healthy_samples[
  hd_tree$order
]

mtnbc_samples <- mtnbc_samples[
  mtnbc_tree$order
]


column_order <- c(
  healthy_samples,
  mtnbc_samples
)

heatmap_mat <- heatmap_mat[
  ,
  column_order,
  drop = FALSE
]


# Convert each gene to a row Z-score
heatmap_z <- t(
  scale(
    t(heatmap_mat)
  )
)

heatmap_z[
  is.na(heatmap_z)
] <- 0


# Use gene symbols on the right side
gene_labels <- sig_df$gene_symbol

missing_symbol <- is.na(gene_labels) |
  gene_labels == ""

gene_labels[missing_symbol] <-
  sig_df$ensembl_gene_id[missing_symbol]

rownames(heatmap_z) <- make.unique(
  gene_labels
)


# Clustering

# Row clustering
row_tree <- hclust(
  dist(heatmap_z),
  method = "average"
)


# Global column clustering
#
# This is used to draw the top dendrogram. If the two biological
# groups separate naturally, it should resemble the paper's tree
# reasonably closely.
column_tree <- hclust(
  dist(
    t(heatmap_z)
  ),
  method = "average"
)


# Bottom HD / mTNBC labels

n_hd <- length(
  healthy_samples
)

n_mtnbc <- length(
  mtnbc_samples
)


bottom_group_annotation <- HeatmapAnnotation(
  
  group = AnnotationFunction(
    
    fun = function(index) {
      
      n <- length(index)
      
      # Position separating the two groups
      hd_end <- n_hd / n
      
      # Horizontal bar under HDs
      grid.lines(
        x = unit(
          c(
            0,
            hd_end - 0.01
          ),
          "npc"
        ),
        y = unit(
          c(
            0.75,
            0.75
          ),
          "npc"
        ),
        gp = gpar(
          col = "black",
          lwd = 2
        )
      )
      
      # Small vertical end marks
      grid.lines(
        x = unit(
          c(0, 0),
          "npc"
        ),
        y = unit(
          c(
            0.60,
            0.90
          ),
          "npc"
        )
      )
      
      grid.lines(
        x = unit(
          c(
            hd_end - 0.01,
            hd_end - 0.01
          ),
          "npc"
        ),
        y = unit(
          c(
            0.60,
            0.90
          ),
          "npc"
        )
      )
      
      # Horizontal bar under mTNBC
      grid.lines(
        x = unit(
          c(
            hd_end + 0.01,
            1
          ),
          "npc"
        ),
        y = unit(
          c(
            0.75,
            0.75
          ),
          "npc"
        ),
        gp = gpar(
          col = "black",
          lwd = 2
        )
      )
      
      grid.lines(
        x = unit(
          c(
            hd_end + 0.01,
            hd_end + 0.01
          ),
          "npc"
        ),
        y = unit(
          c(
            0.60,
            0.90
          ),
          "npc"
        )
      )
      
      grid.lines(
        x = unit(
          c(1, 1),
          "npc"
        ),
        y = unit(
          c(
            0.60,
            0.90
          ),
          "npc"
        )
      )
      
      
      # Group names
      grid.text(
        "HDs",
        x = unit(
          hd_end / 2,
          "npc"
        ),
        y = unit(
          0.15,
          "npc"
        ),
        gp = gpar(
          fontsize = 14
        )
      )
      
      grid.text(
        "mTNBC",
        x = unit(
          hd_end +
            (1 - hd_end) / 2,
          "npc"
        ),
        y = unit(
          0.15,
          "npc"
        ),
        gp = gpar(
          fontsize = 14
        )
      )
    },
    
    height = unit(
      9,
      "mm"
    )
  ),
  
  show_annotation_name = FALSE
)


# Color scale
heatmap_colors <- colorRamp2(
  c(
    -2,
    0,
    2
  ),
  c(
    "royalblue3",
    "white",
    "red"
  )
)


# Build heatmap
ht <- Heatmap(
  
  heatmap_z,
  
  name = "Z-score",
  
  col = heatmap_colors,
  
  cluster_rows = row_tree,
  
  # Use the global column dendrogram
  cluster_columns = column_tree,
  
  show_row_names = TRUE,
  show_column_names = FALSE,
  
  row_names_side = "right",
  
  row_names_gp = gpar(
    fontsize = 4
  ),
  
  bottom_annotation =
    bottom_group_annotation,
  
  column_title =
    "Neutrophil transcriptomics",
  
  column_title_gp = gpar(
    fontsize = 16,
    fontface = "bold"
  ),
  
  heatmap_legend_param = list(
    title = NULL,
    at = c(
      -2,
      -1,
      0,
      1,
      2
    ),
    labels = c(
      "-2",
      "-1",
      "0",
      "1",
      "2"
    )
  ),
  
  row_dend_width = unit(
    25,
    "mm"
  ),
  
  column_dend_height = unit(
    13,
    "mm"
  ),
  
  border = FALSE
)

# Save PNG

png(
  file.path(
    "results",
    "Figure5a_heatmap.png"
  ),
  width = 1300,
  height = 1700,
  res = 200
)

draw(
  ht,
  heatmap_legend_side = "left",
  padding = unit(
    c(
      10,
      10,
      10,
      20
    ),
    "mm"
  )
)

dev.off()


# Save PDF

pdf(
  file.path(
    "results",
    "Figure5a_heatmap.pdf"
  ),
  width = 5.5,
  height = 8.5
)

draw(
  ht,
  heatmap_legend_side = "left",
  padding = unit(
    c(
      10,
      10,
      10,
      20
    ),
    "mm"
  )
)

dev.off()


# Figure 5b - Volcano plot


volcano_df <- res_df[
  !is.na(res_df$pvalue),
  ,
  drop = FALSE
]


# DESeq2 p-value transformed for the volcano y-axis
volcano_df$neg_log10_p <- -log10(
  pmax(
    volcano_df$pvalue,
    .Machine$double.xmin
  )
)



# These cutoffs correspond to the visible dashed lines:
#
# p < 0.05
# log2FC > 1   = red
# log2FC < -1  = blue
#
# The formal DEG list remains based on padj < 0.05.
volcano_df$plot_group <- "Other"

volcano_df$plot_group[
  volcano_df$pvalue < 0.05 &
    volcano_df$log2FoldChange > 1
] <- "Up"

volcano_df$plot_group[
  volcano_df$pvalue < 0.05 &
    volcano_df$log2FoldChange < -1
] <- "Down"


# Paper gene labels
# These are annotation labels visible in the published panel.
# They affect only the text shown on the figure and do not
# affect any differential-expression calculations.
paper_genes <- c(
  "SCN9A",
  "CD177",
  "TSPO",
  "ZNF250",
  "MTRR",
  "KL",
  "FCRL1",
  "CLU",
  "TOX",
  "MMP24",
  "SIGLEC17P",
  "GRAMD1C",
  "CCL3L3",
  "CCL4L2"
)


label_df <- volcano_df[
  volcano_df$gene_label %in% paper_genes,
  ,
  drop = FALSE
]

label_df <- label_df[
  !duplicated(
    label_df$gene_label
  ),
]


# Report labels that are not available in the deposited matrix
missing_paper_genes <- setdiff(
  paper_genes,
  label_df$gene_label
)

if (length(missing_paper_genes) > 0) {
  
  cat(
    "\nPaper labels not found in this dataset:\n"
  )
  
  print(
    missing_paper_genes
  )
}



# Qlucore and DESeq2 might not give identical volcano coordinates.
# To keep the visual panel close to the published 0-6 y-axis,
# values above 6 are capped ONLY for plotting.
# Original p-values and -log10(p) values remain unchanged in
# res_df and all saved result tables.
volcano_df$plot_y <- pmin(
  volcano_df$neg_log10_p,
  6
)

label_df$plot_y <- pmin(
  label_df$neg_log10_p,
  6
)


# Draw volcano plot
volcano_plot <- ggplot() +
  
  # Background genes:
  # white center with black outline
  geom_point(
    data = volcano_df[
      volcano_df$plot_group == "Other",
    ],
    aes(
      x = log2FoldChange,
      y = plot_y
    ),
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.35,
    size = 1.8
  ) +
  
  # Downregulated region:
  # blue fill with black outline
  geom_point(
    data = volcano_df[
      volcano_df$plot_group == "Down",
    ],
    aes(
      x = log2FoldChange,
      y = plot_y
    ),
    shape = 21,
    fill = "royalblue4",
    color = "black",
    stroke = 0.3,
    size = 2
  ) +
  
  # Upregulated region:
  # red fill with black outline
  geom_point(
    data = volcano_df[
      volcano_df$plot_group == "Up",
    ],
    aes(
      x = log2FoldChange,
      y = plot_y
    ),
    shape = 21,
    fill = "tomato4",
    color = "black",
    stroke = 0.3,
    size = 2
  ) +
  
  # Fold-change cutoffs
  geom_vline(
    xintercept = c(
      -1,
      1
    ),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  # p = 0.05 cutoff
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  # Labels shown in the published figure.
  # ggrepel adds the small connector line from text to point.
  geom_text_repel(
    data = label_df,
    aes(
      x = log2FoldChange,
      y = plot_y,
      label = gene_label
    ),
    
    size = 2.8,
    
    color = "black",
    
    segment.color = "black",
    segment.size = 0.4,
    
    min.segment.length = 0,
    
    box.padding = 0.5,
    point.padding = 0.1,
    
    max.overlaps = Inf,
    
    seed = 123
  ) +
  
  # X-axis spacing
  scale_x_continuous(
    breaks = seq(
      -3,
      5,
      by = 1
    )
  ) +
  
  # Paper-style Y-axis
  scale_y_continuous(
    breaks = seq(
      0,
      6,
      by = 0.5
    ),
    limits = c(
      0,
      6.1
    ),
    expand = expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +
  
  coord_cartesian(
    xlim = c(
      -3.5,
      5.5
    )
  ) +
  
  labs(
    title = "Differential gene expression neutrophils",
    subtitle = "mTNBC vs Healthy Donors",
    x = expression(
      "Difference (Log"[2] *
        " Fold Change) of Group Means"
    ),
    y = expression(
      -Log[10](p)
    )
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 14,
      margin = margin(
        b = 15
      )
    ),
    
    axis.title.x = element_text(
      size = 14
    ),
    
    axis.title.y = element_text(
      size = 14
    ),
    
    axis.text = element_text(
      size = 10
    ),
    
    axis.line = element_line(
      linewidth = 0.6,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.5,
      color = "black"
    ),
    
    legend.position = "none"
  )


print(
  volcano_plot
)


# Save volcano as PNG
ggsave(
  filename = file.path(
    "results",
    "Figure5b_volcano.png"
  ),
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# Save volcano as PDF
ggsave(
  filename = file.path(
    "results",
    "Figure5b_volcano.pdf"
  ),
  plot = volcano_plot,
  width = 8,
  height = 6
)


# 7. Record the software environment

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    "results",
    "sessionInfo.txt"
  )
)


# Finished
cat(
  "\nAnalysis complete.\n",
  "Results saved in:\n",
  normalizePath("results"),
  "\n\nFiles created:\n"
)

print(
  list.files("results")
)