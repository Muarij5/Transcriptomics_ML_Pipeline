# ============================================================
# INSTALL ALL REQUIRED PACKAGES — Run once only
# ============================================================
install.packages(c("ggplot2", "ggrepel", "dplyr", "patchwork"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "sva"))

library(ggVennDiagram)
library(GEOquery)
library(DESeq2)
library(sva)
library(ggplot2)
library(dplyr)
library(lme4)
library(httr)
library(readr)
library(org.Hs.eg.db)
library(ggrepel)
library(patchwork)

setwd("D:/Epilepsy_merged")

# ============================================================
# LOAD DATA
# ============================================================
expr_raw     <- read.csv("final_clean_matrix.csv", check.names = FALSE)
gene_symbols <- expr_raw[, 1]
expr_raw     <- expr_raw[, -1]
rownames(expr_raw) <- gene_symbols

metadata <- read.csv("merged_metadata.csv")
rownames(metadata) <- metadata$sample_id

# ============================================================
# ALIGN SAMPLES
# ============================================================
common_samples <- intersect(colnames(expr_raw), metadata$sample_id)
expr_raw       <- expr_raw[, common_samples]
metadata       <- metadata[common_samples, ]

stopifnot(all(colnames(expr_raw) == rownames(metadata)))
cat("Samples aligned:", ncol(expr_raw), "\n")

# ============================================================
# REMOVE OUTLIER TISSUES
# ============================================================
tissues_to_remove <- c("Cortex", "Surgical_Brain", "Brain")

metadata_clean <- metadata[!metadata$tissue %in% tissues_to_remove, ]
expr_raw_clean <- expr_raw[, metadata_clean$sample_id]

cat("Samples after tissue removal:", ncol(expr_raw_clean), "\n")
print(table(metadata_clean$tissue))

# ============================================================
# BROAD CONDITION GROUPING
# ============================================================
metadata_clean$condition_broad <- metadata_clean$condition
metadata_clean$condition_broad[metadata_clean$condition %in%
                                 c("mTLE", "mTLE_HS", "mTLE_nonHS")] <- "mTLE"
metadata_clean$condition_broad[metadata_clean$condition %in%
                                 c("FCD_IIa", "FCD_IIb")]             <- "FCD"
metadata_clean$condition_broad[metadata_clean$condition == "TSC"]     <- "TSC"
metadata_clean$condition_broad[metadata_clean$condition == "Control"] <- "Control"

metadata_clean$condition_broad <- factor(
  metadata_clean$condition_broad,
  levels = c("Control", "mTLE", "FCD", "TSC")
)

cat("\n=== CONDITION_BROAD TABLE ===\n")
print(table(metadata_clean$condition_broad))
# ============================================================
# RE-ALIGN AFTER TISSUE REMOVAL
# ============================================================
common_samples_clean <- intersect(colnames(expr_raw_clean), rownames(metadata_clean))
expr_raw_clean       <- expr_raw_clean[, common_samples_clean]
metadata_clean       <- metadata_clean[common_samples_clean, ]

stopifnot("Mismatch after tissue removal!" =
            identical(colnames(expr_raw_clean), rownames(metadata_clean)))

cat("\n=== SAMPLE ALIGNMENT AFTER TISSUE REMOVAL ===\n")
cat("expr_raw_clean dimensions:  ", nrow(expr_raw_clean), "genes x", ncol(expr_raw_clean), "samples\n")
cat("metadata_clean dimensions:  ", nrow(metadata_clean), "samples x", ncol(metadata_clean), "variables\n")
cat("Column/row name check:       PASSED\n")

# ============================================================
# ETA SQUARED HELPER FUNCTION
# ============================================================
compute_eta_sq <- function(expr_mat, group_labels) {
  sapply(1:ncol(expr_mat), function(i) {
    gene_expr <- expr_mat[, i]
    ss_total  <- sum((gene_expr - mean(gene_expr))^2)
    ss_batch  <- sum(tapply(gene_expr, group_labels, function(x)
      length(x) * (mean(x) - mean(gene_expr))^2))
    return(ss_batch / ss_total)
  })
}

# ============================================================
# VST — visualization only, before batch correction
# ============================================================
dds_raw    <- DESeqDataSetFromMatrix(
  countData = expr_raw_clean,
  colData   = metadata_clean,
  design    = ~ condition_broad
)
vst_mat    <- vst(dds_raw, blind = TRUE)
vst_counts <- assay(vst_mat)

# ============================================================
# PCA BEFORE BATCH CORRECTION
# ============================================================
pca_res       <- prcomp(t(vst_counts), scale. = TRUE)
pca_df        <- as.data.frame(pca_res$x[, 1:2])
var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 2)

pca_df$dataset_id      <- metadata_clean[rownames(pca_df), "dataset_id"]
pca_df$condition_broad <- metadata_clean[rownames(pca_df), "condition_broad"]
pca_df$tissue          <- metadata_clean[rownames(pca_df), "tissue"]
pca_df$platform        <- metadata_clean[rownames(pca_df), "platform"]

p1 <- ggplot(pca_df, aes(PC1, PC2, color = dataset_id)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA Before Batch Correction — By Dataset",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_bw()

p2 <- ggplot(pca_df, aes(PC1, PC2, color = condition_broad)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA Before Batch Correction — By Condition",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_bw()

p3 <- ggplot(pca_df, aes(PC1, PC2, color = tissue)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA Before Batch Correction — By Tissue",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_bw()

p4 <- ggplot(pca_df, aes(PC1, PC2, color = platform)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA Before Batch Correction — By Platform",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_bw()

print(p1); print(p2); print(p3); print(p4)
ggsave("PCA_before_batch_by_dataset.png",   p1, width = 10, height = 7, dpi = 300)
ggsave("PCA_before_batch_by_condition.png", p2, width = 10, height = 7, dpi = 300)
ggsave("PCA_before_batch_by_tissue.png",    p3, width = 10, height = 7, dpi = 300)
ggsave("PCA_before_batch_by_platform.png",  p4, width = 10, height = 7, dpi = 300)

# ============================================================
# ETA² BEFORE BATCH CORRECTION
# ============================================================
gene_var  <- apply(vst_counts, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE)[1:1000])
expr_top  <- t(vst_counts[top_genes, ])

eta_dataset  <- compute_eta_sq(expr_top, metadata_clean$dataset_id)
eta_tissue   <- compute_eta_sq(expr_top, metadata_clean$tissue)
eta_platform <- compute_eta_sq(expr_top, metadata_clean$platform)

cat("\n=== ETA SQUARED BEFORE CORRECTION (VST) ===\n")
cat("dataset_id | Mean:", round(mean(eta_dataset),  4),
    "| Median:", round(median(eta_dataset),  4), "\n")
cat("tissue     | Mean:", round(mean(eta_tissue),   4),
    "| Median:", round(median(eta_tissue),   4), "\n")
cat("platform   | Mean:", round(mean(eta_platform), 4),
    "| Median:", round(median(eta_platform), 4), "\n")

# ComBat → input normalized → output normalized → PCA directly → limma for DEGs
# ComBat-seq → input raw counts → output raw counts → VST → PCA → DESeq2 for DEGs

# ============================================================
# ComBat_seq ON RAW INTEGER COUNTS
# ============================================================
expr_counts <- round(as.matrix(expr_raw_clean))
storage.mode(expr_counts) <- "integer"

cat("\nRunning ComBat_seq on raw counts...\n")
expr_corrected_counts <- ComBat_seq(
  counts   = expr_counts,
  batch    = metadata_clean$dataset_id,
  group    = metadata_clean$condition_broad,
  full_mod = TRUE
)
cat("ComBat_seq complete. Dims:", dim(expr_corrected_counts), "\n")

# ============================================================
# VST ON CORRECTED COUNTS — for PCA, Eta², and ML
# ============================================================
dds_corrected <- DESeqDataSetFromMatrix(
  countData = expr_corrected_counts,
  colData   = metadata_clean,
  design    = ~ condition_broad
)
vst_corrected     <- vst(dds_corrected, blind = FALSE)
vst_corrected_mat <- assay(vst_corrected)

# ============================================================
# ETA² AFTER ComBat_seq
# ============================================================
gene_var2  <- apply(vst_corrected_mat, 1, var)
top_genes2 <- names(sort(gene_var2, decreasing = TRUE)[1:1000])
expr_top2  <- t(vst_corrected_mat[top_genes2, ])

eta_dataset_cb  <- compute_eta_sq(expr_top2, metadata_clean$dataset_id)
eta_tissue_cb   <- compute_eta_sq(expr_top2, metadata_clean$tissue)
eta_platform_cb <- compute_eta_sq(expr_top2, metadata_clean$platform)

cat("\n=== ETA SQUARED AFTER ComBat_seq (VST) ===\n")
cat("dataset_id | Mean:", round(mean(eta_dataset_cb),  4),
    "| Median:", round(median(eta_dataset_cb),  4), "\n")
cat("tissue     | Mean:", round(mean(eta_tissue_cb),   4),
    "| Median:", round(median(eta_tissue_cb),   4), "\n")
cat("platform   | Mean:", round(mean(eta_platform_cb), 4),
    "| Median:", round(median(eta_platform_cb), 4), "\n")

cat("\n=== COMPARISON: Before vs After ComBat_seq ===\n")
cat("           | Before  | After\n")
cat("dataset_id |", round(mean(eta_dataset),  4), "|",
    round(mean(eta_dataset_cb),  4), "\n")
cat("tissue     |", round(mean(eta_tissue),   4), "|",
    round(mean(eta_tissue_cb),   4), "\n")
cat("platform   |", round(mean(eta_platform), 4), "|",
    round(mean(eta_platform_cb), 4), "\n")

# ============================================================
# PCA AFTER ComBat_seq
# ============================================================
zero_var               <- apply(vst_corrected_mat, 1, var) == 0
vst_corrected_filtered <- vst_corrected_mat[!zero_var, ]

pca_cb    <- prcomp(t(vst_corrected_filtered), scale. = TRUE)
pca_cb_df <- as.data.frame(pca_cb$x[, 1:2])
var_cb    <- round(100 * pca_cb$sdev^2 / sum(pca_cb$sdev^2), 2)

pca_cb_df$dataset_id      <- metadata_clean[rownames(pca_cb_df), "dataset_id"]
pca_cb_df$condition_broad <- metadata_clean[rownames(pca_cb_df), "condition_broad"]
pca_cb_df$tissue          <- metadata_clean[rownames(pca_cb_df), "tissue"]
pca_cb_df$platform        <- metadata_clean[rownames(pca_cb_df), "platform"]

p_batch <- ggplot(pca_cb_df, aes(PC1, PC2, color = dataset_id)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "PCA After ComBat_seq — By Dataset",
       x = paste0("PC1 (", var_cb[1], "%)"),
       y = paste0("PC2 (", var_cb[2], "%)")) +
  theme_bw()

p_cond <- ggplot(pca_cb_df, aes(PC1, PC2, color = condition_broad)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "PCA After ComBat_seq — By Condition",
       x = paste0("PC1 (", var_cb[1], "%)"),
       y = paste0("PC2 (", var_cb[2], "%)")) +
  theme_bw()

p_tissue <- ggplot(pca_cb_df, aes(PC1, PC2, color = tissue)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "PCA After ComBat_seq — By Tissue",
       x = paste0("PC1 (", var_cb[1], "%)"),
       y = paste0("PC2 (", var_cb[2], "%)")) +
  theme_bw()

p_platform <- ggplot(pca_cb_df, aes(PC1, PC2, color = platform)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "PCA After ComBat_seq — By Platform",
       x = paste0("PC1 (", var_cb[1], "%)"),
       y = paste0("PC2 (", var_cb[2], "%)")) +
  theme_bw()

print(p_batch); print(p_cond); print(p_tissue); print(p_platform)
ggsave("PCA_after_ComBatseq_by_dataset.png",   p_batch,    width = 10, height = 7, dpi = 300)
ggsave("PCA_after_ComBatseq_by_condition.png", p_cond,     width = 10, height = 7, dpi = 300)
ggsave("PCA_after_ComBatseq_by_tissue.png",    p_tissue,   width = 10, height = 7, dpi = 300)
ggsave("PCA_after_ComBatseq_by_platform.png",  p_platform, width = 10, height = 7, dpi = 300)

# ============================================================
# SAVE CORRECTED COUNTS AND METADATA
# ============================================================
write.csv(as.data.frame(expr_corrected_counts),
          "expr_combatseq_corrected.csv", row.names = TRUE)
write.csv(metadata_clean,
          "metadata_clean.csv", row.names = TRUE)
cat("Saved: expr_combatseq_corrected.csv\n")
cat("Saved: metadata_clean.csv\n")

# ============================================================
# SAVE VST MATRIX FOR ML — samples x genes + condition label
# ============================================================
ml_matrix           <- as.data.frame(t(vst_corrected_mat))
ml_matrix$condition <- metadata_clean[rownames(ml_matrix), "condition_broad"]
write.csv(ml_matrix, "ML_matrix_VST.csv", row.names = TRUE)
cat("Saved: ML_matrix_VST.csv\n")

print(table(metadata_clean$dataset_id, metadata_clean$condition_broad))

# ============================================================
# DESeq2 — filter and run on dds_corrected
# ============================================================
cat("\nFiltering low count genes...\n")
keep          <- rowSums(counts(dds_corrected)) >= 10
dds_corrected <- dds_corrected[keep, ]
cat("Genes after filtering:", nrow(dds_corrected), "\n")

cat("Running DESeq2...\n")
dds_corrected <- DESeq(dds_corrected)
cat("DESeq2 complete. Dims:", dim(dds_corrected), "\n")

# ============================================================
# EXTRACT DEG RESULTS
# ============================================================
res_mTLE_vs_Control <- results(dds_corrected, contrast = c("condition_broad", "mTLE", "Control"))
res_FCD_vs_Control  <- results(dds_corrected, contrast = c("condition_broad", "FCD",  "Control"))
res_TSC_vs_Control  <- results(dds_corrected, contrast = c("condition_broad", "TSC",  "Control"))

cat("\n=== SUMMARY: mTLE vs Control ===\n"); summary(res_mTLE_vs_Control)
cat("\n=== SUMMARY: FCD vs Control ===\n");  summary(res_FCD_vs_Control)
cat("\n=== SUMMARY: TSC vs Control ===\n");  summary(res_TSC_vs_Control)

# ============================================================
# ABLATION STUDY
# ============================================================
deg_list <- list(
  "mTLE_vs_Control" = res_mTLE_vs_Control,
  "FCD_vs_Control"  = res_FCD_vs_Control,
  "TSC_vs_Control"  = res_TSC_vs_Control
)

log2fc_cutoffs <- c(1, 1.5, 2)
padj_cutoffs   <- c(0.05, 0.01, 0.001)

ablation_results <- expand.grid(
  File   = names(deg_list),
  log2FC = c(log2fc_cutoffs, NA),
  padj   = c(padj_cutoffs,   NA)
) %>%
  mutate(Type = "", Total_DEGs = 0, Upregulated = 0, Downregulated = 0)

for(i in 1:nrow(ablation_results)) {
  file_name <- ablation_results$File[i]
  fc        <- ablation_results$log2FC[i]
  p         <- ablation_results$padj[i]
  res       <- as.data.frame(deg_list[[file_name]])
  degs      <- res %>% filter(!is.na(padj))
  
  if(!is.na(fc) & !is.na(p)) {
    degs <- degs %>% filter(abs(log2FoldChange) >= fc & padj <= p)
    ablation_results$Type[i] <- "both"
  } else if(is.na(fc) & !is.na(p)) {
    degs <- degs %>% filter(padj <= p)
    ablation_results$Type[i] <- "padj-only"
  } else if(!is.na(fc) & is.na(p)) {
    degs <- degs %>% filter(abs(log2FoldChange) >= fc)
    ablation_results$Type[i] <- "log2FC-only"
  }
  
  ablation_results$Total_DEGs[i]    <- nrow(degs)
  ablation_results$Upregulated[i]   <- sum(degs$log2FoldChange > 0)
  ablation_results$Downregulated[i] <- sum(degs$log2FoldChange < 0)
}

print(ablation_results)
write.csv(ablation_results, "DEG_Ablation_Study.csv", row.names = FALSE)

# ============================================================
# FINAL DEG FILTERING (padj < 0.05 & |log2FC| > 1)
# ============================================================
DEGs_mTLE_vs_Control <- subset(res_mTLE_vs_Control, padj < 0.05 & abs(log2FoldChange) > 1)
DEGs_FCD_vs_Control  <- subset(res_FCD_vs_Control,  padj < 0.05 & abs(log2FoldChange) > 1)
DEGs_TSC_vs_Control  <- subset(res_TSC_vs_Control,  padj < 0.05 & abs(log2FoldChange) > 1)

cat("=== DEG Counts (padj < 0.05 & |log2FC| > 1) ===\n")
cat("mTLE vs Control:", nrow(DEGs_mTLE_vs_Control), "\n")
cat("FCD  vs Control:", nrow(DEGs_FCD_vs_Control),  "\n")
cat("TSC  vs Control:", nrow(DEGs_TSC_vs_Control),  "\n")

DEGs_mTLE_vs_Control$Direction <- ifelse(DEGs_mTLE_vs_Control$log2FoldChange > 0, "Up", "Down")
DEGs_FCD_vs_Control$Direction  <- ifelse(DEGs_FCD_vs_Control$log2FoldChange  > 0, "Up", "Down")
DEGs_TSC_vs_Control$Direction  <- ifelse(DEGs_TSC_vs_Control$log2FoldChange  > 0, "Up", "Down")
# After your DEG filtering lines
DEGs_mTLE_vs_Control <- as.data.frame(DEGs_mTLE_vs_Control)
DEGs_FCD_vs_Control  <- as.data.frame(DEGs_FCD_vs_Control)
DEGs_TSC_vs_Control  <- as.data.frame(DEGs_TSC_vs_Control)
write.csv(as.data.frame(DEGs_mTLE_vs_Control), "DEGs_mTLE_vs_Control_final.csv", row.names = TRUE)
write.csv(as.data.frame(DEGs_FCD_vs_Control),  "DEGs_FCD_vs_Control_final.csv",  row.names = TRUE)
write.csv(as.data.frame(DEGs_TSC_vs_Control),  "DEGs_TSC_vs_Control_final.csv",  row.names = TRUE)

dim(DEGs_mTLE_vs_Control)
dim(DEGs_FCD_vs_Control)
dim(DEGs_TSC_vs_Control)
cat("\nAll files saved successfully.\n")

# ============================================================
# GLOBAL THEME 
# ============================================================
nature_theme <- theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.65),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(3, "mm"),
    axis.title = element_text(size = 12.5, face = "bold", color = "black"),
    axis.text = element_text(size = 10.5, color = "black"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 13.5, face = "bold", hjust = 0.5,
                              color = "black", margin = margin(b = 12)),
    legend.title = element_blank(),
    legend.text = element_text(size = 9.5),
    legend.key.size = unit(5, "mm"),
    legend.background = element_rect(fill = NA, color = NA),
    plot.margin = margin(12, 15, 12, 15, "mm")
  )

# Colors 
col_up   <- "#C0392B"
col_down <- "#2471A3"
col_ns   <- "#7F8C8E"

# ============================================================
# HELPER: Classify DEGs
# ============================================================
annotate_degs <- function(res_obj, padj_cut = 0.05, lfc_cut = 1) {
  as.data.frame(res_obj) %>%
    tibble::rownames_to_column("gene") %>%
    filter(!is.na(padj), !is.na(log2FoldChange), !is.na(baseMean)) %>%
    mutate(
      Status = case_when(
        padj < padj_cut & log2FoldChange > lfc_cut  ~ "Up",
        padj < padj_cut & log2FoldChange < -lfc_cut ~ "Down",
        TRUE ~ "NS"
      ),
      Status = factor(Status, levels = c("Up", "Down", "NS"))
    )
}

# ============================================================
# VOLCANO PLOT FUNCTION (Improved)
# ============================================================
make_volcano <- function(res_obj, title, padj_cut = 0.05, lfc_cut = 1) {
  df <- annotate_degs(res_obj, padj_cut, lfc_cut) %>%
    mutate(neglog10p = pmin(-log10(padj), 250))   # cap extremes
  
  n_up  <- sum(df$Status == "Up")
  n_down <- sum(df$Status == "Down")
  n_ns  <- sum(df$Status == "NS")
  
  max_lfc <- ceiling(max(abs(df$log2FoldChange), na.rm = TRUE) + 0.6)
  
  ggplot(df, aes(x = log2FoldChange, y = neglog10p)) +
    # Threshold lines (subtle)
    geom_hline(yintercept = -log10(padj_cut),
               linetype = "dashed", color = "grey55", linewidth = 0.4) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut),
               linetype = "dashed", color = "grey55", linewidth = 0.4) +
    
    # NS points first (background)
    geom_point(data = filter(df, Status == "NS"),
               color = col_ns, size = 0.45, alpha = 0.35) +
    
    # Significant points on top
    geom_point(data = filter(df, Status == "Down"),
               color = col_down, size = 0.95, alpha = 0.88) +
    geom_point(data = filter(df, Status == "Up"),
               color = col_up, size = 0.95, alpha = 0.88) +
    
    scale_x_continuous(
      limits = c(-max_lfc, max_lfc),
      expand = expansion(mult = 0.03),
      breaks = scales::pretty_breaks(n = 7)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05)),
      breaks = scales::pretty_breaks(n = 6)
    ) +
    
    # Clean internal legend (top-left)
    annotate("point", x = -max_lfc + 0.8, y = max(df$neglog10p) * 0.96,
             color = col_up, size = 2.4) +
    annotate("text", x = -max_lfc + 1.4, y = max(df$neglog10p) * 0.96,
             label = paste0("Up (", n_up, ")"), hjust = 0, size = 3.4) +
    
    annotate("point", x = -max_lfc + 0.8, y = max(df$neglog10p) * 0.88,
             color = col_down, size = 2.4) +
    annotate("text", x = -max_lfc + 1.4, y = max(df$neglog10p) * 0.88,
             label = paste0("Down (", n_down, ")"), hjust = 0, size = 3.4) +
    
    annotate("point", x = -max_lfc + 0.8, y = max(df$neglog10p) * 0.80,
             color = col_ns, size = 2.4) +
    annotate("text", x = -max_lfc + 1.4, y = max(df$neglog10p) * 0.80,
             label = paste0("NS (", n_ns, ")"), hjust = 0, size = 3.4) +
    
    labs(
      title = title,
      x = expression(log[2] ~ "Fold Change"),
      y = expression(-log[10] ~ "(adjusted" ~ italic(p) ~ "value)")
    ) +
    nature_theme +
    theme(legend.position = "none")
}

# ============================================================
# MA PLOT FUNCTION (Improved)
# ============================================================
make_MA <- function(res_obj, title, padj_cut = 0.05, lfc_cut = 1) {
  df <- annotate_degs(res_obj, padj_cut, lfc_cut) %>%
    mutate(log10mean = log10(baseMean + 1))
  
  n_up  <- sum(df$Status == "Up")
  n_down <- sum(df$Status == "Down")
  n_ns  <- sum(df$Status == "NS")
  
  max_lfc <- ceiling(max(abs(df$log2FoldChange), na.rm = TRUE) + 0.4)
  max_x   <- max(df$log10mean, na.rm = TRUE) * 1.02
  
  ggplot(df, aes(x = log10mean, y = log2FoldChange)) +
    # Reference lines
    geom_hline(yintercept = 0, color = "black", linewidth = 0.55) +
    geom_hline(yintercept = c(-lfc_cut, lfc_cut),
               linetype = "dashed", color = "grey55", linewidth = 0.4) +
    
    # NS points
    geom_point(data = filter(df, Status == "NS"),
               color = col_ns, size = 0.45, alpha = 0.35) +
    
    # Significant points
    geom_point(data = filter(df, Status == "Down"),
               color = col_down, size = 0.95, alpha = 0.88) +
    geom_point(data = filter(df, Status == "Up"),
               color = col_up, size = 0.95, alpha = 0.88) +
    
    scale_x_continuous(
      expand = expansion(mult = 0.03),
      breaks = scales::pretty_breaks(n = 6)
    ) +
    scale_y_continuous(
      limits = c(-max_lfc, max_lfc),
      expand = expansion(mult = 0.03),
      breaks = scales::pretty_breaks(n = 7)
    ) +
    
    # Clean internal legend (top-right)
    annotate("point", x = max_x * 0.68, y = max_lfc * 0.93,
             color = col_up, size = 2.4) +
    annotate("text", x = max_x * 0.72, y = max_lfc * 0.93,
             label = paste0("Up (", n_up, ")"), hjust = 0, size = 3.4) +
    
    annotate("point", x = max_x * 0.68, y = max_lfc * 0.80,
             color = col_down, size = 2.4) +
    annotate("text", x = max_x * 0.72, y = max_lfc * 0.80,
             label = paste0("Down (", n_down, ")"), hjust = 0, size = 3.4) +
    
    annotate("point", x = max_x * 0.68, y = max_lfc * 0.67,
             color = col_ns, size = 2.4) +
    annotate("text", x = max_x * 0.72, y = max_lfc * 0.67,
             label = paste0("NS (", n_ns, ")"), hjust = 0, size = 3.4) +
    
    labs(
      title = title,
      x = expression(log[10] ~ "Mean Expression"),
      y = expression(log[2] ~ "Fold Change")
    ) +
    nature_theme +
    theme(legend.position = "none")
}

# ============================================================
# SAVE HELPER 
# ============================================================
save_plot <- function(p, name, w = 120, h = 110) {
  ggsave(paste0(name, ".pdf"), plot = p,
         width = w, height = h, units = "mm",
         dpi = 300, device = cairo_pdf)
  ggsave(paste0(name, ".png"), plot = p,
         width = w, height = h, units = "mm",
         dpi = 300, bg = "white")
  cat("Saved:", name, "\n")
}

# ============================================================
# BUILD PLOTS
# ============================================================
v1 <- make_volcano(res_mTLE_vs_Control, "mTLE vs. Control")
v2 <- make_volcano(res_FCD_vs_Control, "FCD vs. Control")
v3 <- make_volcano(res_TSC_vs_Control, "TSC vs. Control")

m1 <- make_MA(res_mTLE_vs_Control, "mTLE vs. Control")
m2 <- make_MA(res_FCD_vs_Control, "FCD vs. Control")
m3 <- make_MA(res_TSC_vs_Control, "TSC vs. Control")

# ============================================================
# SAVE INDIVIDUAL PLOTS
# ============================================================
save_plot(v1, "Volcano_mTLE_vs_Control", w = 120, h = 110)
save_plot(v2, "Volcano_FCD_vs_Control", w = 120, h = 110)
save_plot(v3, "Volcano_TSC_vs_Control", w = 120, h = 110)

save_plot(m1, "MA_mTLE_vs_Control", w = 120, h = 110)
save_plot(m2, "MA_FCD_vs_Control", w = 120, h = 110)
save_plot(m3, "MA_TSC_vs_Control", w = 120, h = 110)

# ============================================================
# COMBINED PANEL 
# ============================================================
combined <- (v1 | v2 | v3) / (m1 | m2 | m3) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 14, face = "bold", family = "Helvetica"))
  )

save_plot(combined, "Figure_DEG_Panel", w = 183, h = 200)

cat("\n=== All plots updated and saved. ===\n")
# ============================================================
# 1. EXTRACT GENE LISTS
# ============================================================

genes_mTLE <- rownames(DEGs_mTLE_vs_Control)
genes_FCD  <- rownames(DEGs_FCD_vs_Control)
genes_TSC  <- rownames(DEGs_TSC_vs_Control)

gene_sets <- list(
  "mTLE vs Control" = genes_mTLE,
  "FCD vs Control"  = genes_FCD,
  "TSC vs Control"  = genes_TSC
)

# ============================================================
# 2. COMMON GENES ACROSS ALL THREE CONDITIONS
# ============================================================

common_genes <- Reduce(intersect, gene_sets)
cat("=== Common DEGs across ALL 3 conditions ===\n")
cat("Total common genes:", length(common_genes), "\n")
print(common_genes)

# Save common genes (was referencing non-existent 'common_genes_filtered')
write.csv(
  data.frame(Gene = common_genes),
  "Common_DEGs_all3_conditions.csv",
  row.names = FALSE
)
cat("Saved: Common_DEGs_all3_conditions.csv\n")

# ============================================================
# 3. PAIRWISE OVERLAPS (for reporting)
# ============================================================

common_mTLE_FCD <- intersect(genes_mTLE, genes_FCD)
common_mTLE_TSC <- intersect(genes_mTLE, genes_TSC)
common_FCD_TSC  <- intersect(genes_FCD,  genes_TSC)

cat("\n=== Pairwise Overlaps ===\n")
cat("mTLE ∩ FCD:", length(common_mTLE_FCD), "\n")
cat("mTLE ∩ TSC:", length(common_mTLE_TSC), "\n")
cat("FCD  ∩ TSC:", length(common_FCD_TSC),  "\n")
cat("mTLE ∩ FCD ∩ TSC:", length(common_genes), "\n")
# ============================================================
# 4.VENN DIAGRAM
# ============================================================

venn_plot <- ggVennDiagram(
  gene_sets,
  label         = "count",
  label_alpha   = 0,
  edge_size     = 1.2,
  set_size      = 6.5
) +
  scale_fill_gradientn(
    colours   = c("#EAF4FB", "#7EC8E3", "#1A6FA8", "#0A3055"),
    name      = "Gene\nCount",
    guide     = guide_colorbar(
      barwidth      = 1.2,
      barheight     = 8,
      title.hjust   = 0.5,
      title.theme   = element_text(size = 11, face = "bold", family = "Arial"),
      label.theme   = element_text(size = 10, family = "Arial")
    )
  ) +
  scale_color_manual(values = c(
    "mTLE vs Control" = "#1A6FA8",
    "FCD vs Control"  = "#E05C2A",
    "TSC vs Control"  = "#2EAA6E"
  )) +
  labs(
    title    = "Differentially Expressed Genes Overlap",
    subtitle = "mTLE, FCD & TSC vs Control  |  padj < 0.05  &  |log₂FC| > 1"
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title    = element_text(
      size   = 16, face   = "bold",
      hjust  = 0.5, margin = margin(b = 6),
      family = "Arial"
    ),
    plot.subtitle = element_text(
      size   = 11, hjust  = 0.5,
      color  = "#444444", margin = margin(b = 14),
      family = "Arial"
    ),
    plot.margin   = margin(20, 30, 20, 30),
    legend.position = "right",
    legend.title  = element_text(size = 11, face = "bold", family = "Arial"),
    legend.text   = element_text(size = 10, family = "Arial"),
    plot.background = element_rect(fill = "white", color = NA)
  )

# ============================================================
# 5. SAVE AT 300 DPI
# ============================================================

ggsave(
  filename  = "Venn_DEGs_publication.png",
  plot      = venn_plot,
  width     = 9,
  height    = 7,
  dpi       = 300,
  units     = "in",
  bg        = "white"
)

cat("\nVenn diagram saved: PNG (300 DPI)")
cat("Common genes CSV saved.\n")


