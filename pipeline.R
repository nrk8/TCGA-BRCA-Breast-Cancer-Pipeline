#install.packages("BiocManager")
BiocManager::install(c("recount3", "SummarizedExperiment"), update = FALSE, ask = FALSE)
library(recount3)
library(SummarizedExperiment)
# Tells R to create the cache folder automatically without asking (yes/no)
options(recount3_cache = "GDCdata")

# This downloads the full breast cancer dataset in one step
se_brca <- create_rse(subset(available_projects(), project == "BRCA" & project_type == "data_sources"))

#1. Properly scale the counts using recount3's transformer
assay(se_brca, "counts") <- compute_read_counts(se_brca)
# Extract correctly scaled expression matrix
gene_expression <- assay(se_brca, "counts")

# 2. Extract medical labels sheet using the correct 'tcga.' prefix
sample_labels <- colData(se_brca)$tcga.gdc_cases.samples.sample_type

# 3. Create a clean list of just the Tumor Tissue columns
tumor_counts <- gene_expression[, sample_labels == "Primary Tumor", drop = FALSE]

# 4. Create a clean list of just the Healthy Tissue columns
healthy_counts <- gene_expression[, sample_labels == "Solid Tissue Normal", drop = FALSE]

print(dim(tumor_counts))
print(dim(healthy_counts))

# 1. Calculate the average activity level for each gene in Tumors and Healthy cells
mean_tumor <- rowMeans(tumor_counts)
mean_healthy <- rowMeans(healthy_counts)

# 2. Calculate the "Fold Change" (how many times more active a gene is in cancer)
# use log2 transformation to make the numbers easy to compare
log2_fold_change <- log2((mean_tumor + 1) / (mean_healthy + 1))

# 3. Combine results into a clean spreadsheet
gene_mapping_results <- data.frame(
  Gene_ID = rownames(tumor_counts),
  Tumor_Average = mean_tumor,
  Healthy_Average = mean_healthy,
  Log2FoldChange = log2_fold_change
)

# 4. Sort the list to put the most active cancer genes at the top
highly_active_in_cancer <- gene_mapping_results[order(-gene_mapping_results$Log2FoldChange), ]

# Print the top 5 most highly active genes in tumor cells
head(highly_active_in_cancer, 5)

# Making a clean table of the active genes in script editor
table_1 <- head(highly_active_in_cancer, 5)
View(table_1)

# 1. Install and load limma for statistical analysis
library(limma)
library(ggplot2)

# 2. Combine tumor and healthy data into one master matrix
combined_matrix <- cbind(tumor_counts, healthy_counts)

# 3. Create a matching label sheet telling R which columns are Tumor vs Healthy
conditions <- c(rep("Tumor", ncol(tumor_counts)), rep("Healthy", ncol(healthy_counts)))
design <- model.matrix(~0 + conditions)
colnames(design) <- c("Healthy", "Tumor")

# 4. Run the linear modeling and empirical Bayes statistics 
fit <- lmFit(combined_matrix, design)
contrast_matrix <- makeContrasts(Tumor - Healthy, levels = design)
fit_contrast <- contrasts.fit(fit, contrast_matrix)
fit_bayes <- eBayes(fit_contrast)

# 5. Extract the final results table with P-values included
final_dge_results <- topTable(fit_bayes, coef = 1, number = Inf)
final_dge_results$Gene_ID <- rownames(final_dge_results)

table_2 <- head(final_dge_results[final_dge_results$logFC > 2 & final_dge_results$adj.P.Val < 0.05, ], 5)
colnames(table_2)[colnames(table_2) == "logFC"] <- "Log2FoldChange"

# 6. Mark genes as significantly Up-regulated or Down-regulated
final_dge_results$Significance <- "Not Significant"
final_dge_results$Significance[final_dge_results$logFC > 2 & final_dge_results$adj.P.Val < 0.05] <- "Highly Active in Cancer"
final_dge_results$Significance[final_dge_results$logFC < -2 & final_dge_results$adj.P.Val < 0.05] <- "Highly Active in Healthy"

# 7. Create a visual Volcano Plot graph
ggplot(final_dge_results, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("red", "blue", "grey")) +
  theme_minimal() +
  labs(title = "RNA-Seq Differential Gene Expression Mapping",
       subtitle = "Breast Cancer Tumors vs. Healthy Tissue Controls",
       x = "Activity Shift Scale (Log2 Fold Change)",
       y = "Statistical Significance Score (-log10 Adjusted P-Value)")


# 1. Choose top 5 most active genes from dataset
graph_data <- table_2
View(table_2)

# 2. Add biological names for text labels
graph_data$Gene_Label <- c("MUC21 (Immune Shield)", 
                           "GRHL2 (Growth Switch)", 
                           "LINC00992 (Cancer Marker)", 
                           "UTY (DNA Packager)", 
                           "LINC01133 (Tumor Supervisor)")

# 3. Draw the labeled bar chart
ggplot(graph_data, aes(x = reorder(Gene_Label, Log2FoldChange), y = Log2FoldChange, fill = Gene_Label)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  coord_flip() + # Flips bars sideways so names are easy to read
  theme_minimal() +
  theme(legend.position = "none", # Removes unneeded color legend box
        plot.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 11)) + 
  scale_fill_brewer(palette = "Reds") + 
  labs(
    title = "Top 5 Most Highly Active Genes in Cancer Cells",
    subtitle = "Calculated by Activity Boost (Log2 Fold Change) vs. Healthy Tissue",
    x = "Specific Biological Gene Name",
    y = "Activity Intensity Boost Scale"
  )

# Install and load the pheatmap library 
if (!requireNamespace("pheatmap")) install.packages("pheatmap")
library(pheatmap)

# 1. Log-Transform data for unbiased clustering 
# add a pseudo-count of 1 to avoid running log2 on absolute zeros
log_matrix <- log2(combined_matrix + 1)


# GRAPH 3: PRINCIPAL COMPONENT ANALYSIS (PCA) PLOT

# Run PCA on the transposed matrix (samples as rows, genes as columns)
# Select the top 500 most variable genes to capture clean biological grouping
var_genes <- order(rowVars(log_matrix), decreasing = TRUE)[1:500]
pca_input <- t(log_matrix[var_genes, ])
pca_results <- prcomp(pca_input, scale. = TRUE)

# Shape the PCA coordinates into a clean plotting data frame
pca_df <- as.data.frame(pca_results$x)
pca_df$Condition <- conditions # Uses existing "Tumor" and "Healthy" labels

# Calculate percentage variance explained by each main axis
percent_var <- round(100 * (pca_results$sdev^2 / sum(pca_results$sdev^2)))

# Draw the PCA scatterplot
ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal() +
  labs(
    title = "Principal Component Analysis (PCA) Cohort Overviews",
    subtitle = "Evaluating Global Variance Patterns Across 1,241 Patient Samples",
    x = paste0("PC1: ", percent_var[1], "% Variance Explained"),
    y = paste0("PC2: ", percent_var[2], "% Variance Explained")
  )


# GRAPH 4: HIERARCHICAL CLUSTERING HEATMAP

# Log-transform data for unbiased clustering
log_matrix <- log2(combined_matrix + 1)

# Select the top 30 significant genes
top_30_hits <- head(final_dge_results[order(final_dge_results$adj.P.Val), ][1:30, ])

# Get corresponding readable gene symbols using manual mapping framework as a temporary bridge to avoid biomart connection lag
heatmap_matrix <- log_matrix[top_30_hits$Gene_ID, , drop = FALSE ]

# Create a readable name vector for the plot rows
readable_labels <- rownames(heatmap_matrix)

# Replace specific known gene IDs with clear names
readable_labels <- gsub("ENSG00000215269", "MUC21", readable_labels)
readable_labels <- gsub("ENSG00000132446", "GRHL2", readable_labels)
readable_labels <- gsub("ENSG00000224902", "LINC00992", readable_labels)

rownames(heatmap_matrix) <- readable_labels

# Setup annotations
annotation_df <- data.frame(Cohort = conditions)
rownames(annotation_df) <- colnames(heatmap_matrix)
annotation_colors <- list(Cohort = c(Healthy = "royalblue3", Tumor = "firebrick2"))

# Draw heatmap with expanded margins
pheatmap(heatmap_matrix,
         scale = "row",
         show_colnames = FALSE,
         show_rownames = TRUE, # shows row gene text clearly
         fontsize_row = 8, # sizes down text slightly so nothing overlaps
         annotation_col = annotation_df,
         annotation_colors = annotation_colors,
         main = "Hierarchical Clustering: Top 30 Variant Driver Transcripts",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         margin = c(5, 10, 5, 15)) # stop text cutting off

