# TCGA-BRCA-Breast-Cancer-Pipeline
An automated RNA-Seq differential gene expression pipeline streaming 1,241 patient samples from the GDC cloud using recount3 and limma.

Methods

Breast cancer RNA-Seq data (TCGA-BRCA) was streamed via recount3 to compare Primary Tumors (n = 1,127) against Normal Controls (n = 114). Raw counts were computed using compute_read_counts() across 63,856 genes. Differential gene expression was executed using the limma package (lmFit and eBayes) to account for count variance. Significant targets were isolated using an absolute threshold of log_2 Fold Change > 2 and an adjusted p-value < 0.05. Top variants were log-scaled and visualized using pheatmap hierarchical clustering.

Results

The limma pipeline isolated a highly significant malignancy signature profile. Hierarchical clustering of the top 30 variant transcripts achieved independent separation between the healthy control and primary tumor cohorts. Three primary oncogenic drivers were visually isolated on the final heatmap: MUC21 and GRHL2 demonstrated uniform, intense over-activation in tumor tissues, while LINC00992 displayed a sharp switch, transitioning from complete silence in healthy controls to aggressive expression in the malignant cohort.

Discussion

The mapped transcripts reveal key mechanisms of tumor progression and clinical utility. MUC21 hyperactivation suggests the creation of a physical sugar shield that blocks immune cell activity. GRHL2 upregulation indicates the activation of a master epithelial-mesenchymal transition (EMT) switch that drives cell detachment and metastasis. Most notably, the binary profile of LINC00992 marks it as an excellent non-invasive biomarker candidate for developing early-stage liquid biopsy screening tests.
