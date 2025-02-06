

# Setting working directory and reading file
setwd("~/bulkcell/Triple-Negative-Breast-Cancer")

metadata <- read.csv("GSE264108_readcounts.txt", sep="\t", row.names = 1)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(pheatmap)

#Extracting gene names from the dataset as external gene id is provided in raw data
gene_names <- metadata$external_gene_id
names(gene_names) <- rownames(metadata)

#Extracting counts and column names to create coldata as there are 14 samples in the dataset
counts <- metadata[,1:14]
col_names <- colnames(counts)

# Extracting condition and naming Tonic(stage I-III) and TNB as mTNBC and Healthy as Healthy
conditions <- ifelse(grepl("Tonic|TNB", col_names), "mTNBC", "Healthy")
print(conditions)

# Creating a colData dataframe based on sample conditions
coldata <- data.frame(
  condition = conditions,
  row.names = col_names
)

# Verifying of all rownames of coldata matches column names of count
all(rownames(coldata) %in% colnames(counts))

# Verifying if rownames of coldata are in same order as column names of count
all(rownames(coldata) == colnames(counts))

#Creating a dds dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,            
  colData = coldata,        
  design = ~ condition        
)

# Plotting the distribution of gene counts
gene_counts <- rowSums(counts(dds))
ggplot(data.frame(gene_counts), aes(x = gene_counts)) +
  geom_histogram(bins = 100, fill = "blue", alpha = 0.7) +
  scale_x_log10() +  # Log-transform x-axis for better visualization
  labs(title = "Distribution of Gene Counts",
       x = "Total Counts per Gene (log10 scale)",
       y = "Frequency") +
  theme_minimal()


# Setting a manual reference as mTNBC as otherwise it will choose automatically
dds$condition <- relevel(dds$condition, ref = "mTNBC")

#Running DESeq2
dds <-  DESeq(dds)

# Comparing mTNBC vs Healthy using reuslts
res_mTNBC_vs_Healthy <- results(dds, contrast = c("condition", "mTNBC", "Healthy"))

# Extracting gene names in the result for the plot later
res_mTNBC_vs_Healthy$gene_name <- gene_names[rownames(res_mTNBC_vs_Healthy)]
summary(res_mTNBC_vs_Healthy)

#MA plot to visualize the upregulated and downregulated genes
plotMA(res_mTNBC_vs_Healthy, ylim=c(-2,2))

# Converting the results to a data frame and removing NA values and adding a significance column
res_df <- as.data.frame(na.omit(res_mTNBC_vs_Healthy))
res_df <- res_df %>%
  mutate(Significance = case_when(
    padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
    padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
    TRUE ~ "Not significant"
  ))

# Identify top 10 upregulated genes
top_upregulated <- res_df %>%
  filter(Significance == "Upregulated") %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)

# Identify top 10 downregulated genes
top_downregulated <- res_df %>%
  filter(Significance == "Downregulated") %>%
  arrange(log2FoldChange) %>%
  head(10)

# Combining top upregulated and downregulated genes for labeling
top_genes <- bind_rows(top_upregulated, top_downregulated)

# Create the Volcano Plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Significance)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "Differential Gene Expression Neutophils
    mTNBC vs Healthy Donors",
    x = "Difference (Log2 Fold Change) of Group Means",
    y = "-Log10 p-value"
  ) +
  geom_text_repel(
    data = top_genes, 
    aes(label = gene_name),
    box.padding = 0.5, 
    max.overlaps = Inf, 
    size = 4, 
    color = "black"  
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), 
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12) 
  )

# Displaying the Volcano plot
print(volcano_plot)



#Set-up for Heatmap Generation--------------------------------------------------

# Filtering for only significant DEGs (upregulated and downregulated)
significant_genes <- res_df %>%
  filter(Significance %in% c("Upregulated", "Downregulated"))

# Extracting normalized counts for all significant DEGs
significant_counts <- normalized_counts[rownames(significant_genes), ]

# Replacing Ensembl IDs with gene names in the row names of the heatmap data
rownames(significant_counts) <- gene_names[rownames(significant_counts)]

# Reordering columns to have "Healthy" first, followed by "mTNBC"
healthy_samples <- colnames(significant_counts)[coldata$condition == "Healthy"]
mtnbc_samples <- colnames(significant_counts)[coldata$condition == "mTNBC"]

# Reordering the counts matrix
significant_counts <- significant_counts[, c(healthy_samples, mtnbc_samples)]

# Reordering the coldata to match the new column order
coldata_reordered <- coldata[c(healthy_samples, mtnbc_samples), , drop = FALSE]

# Creating a column annotation for the heatmap
annotation_col <- data.frame( Condition = coldata_reordered$condition,
                              row.names = colnames(significant_counts))

# Defining colors for the conditions
ann_colors <- list(Condition = c(Healthy = "blue", mTNBC = "red"))

# Creating the heatmap with a color range from -2 to 2
pheatmap(
  mat = significant_counts,
  scale = "row",  
  cluster_rows = TRUE, 
  cluster_cols = FALSE,  
  show_rownames = TRUE, 
  show_colnames = FALSE, 
  annotation_col = annotation_col,  
  annotation_colors = ann_colors,  
  color = colorRampPalette(c("blue", "white", "red"))(100),  
  breaks = seq(-2, 2, length.out = 100), 
  fontsize_row = 8,  
  gaps_col = sum(annotation_col$Condition == "Healthy"), 
  main = "Heatmap of Neutrophil Transcriptomics (Healthy vs mTNBC)"  # Add a title
)



