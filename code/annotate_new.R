############################
# 1️⃣ Load Libraries
############################
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR")

library(biomaRt)
library(data.table)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(edgeR)
options(timeout = 600)

setwd("D:/Epilepsy_merged")

# ==================================================================
# HELPER: Reproducible duplicate resolution by MEDIAN expression
# ==================================================================
dedup_by_median <- function(df, sym_col = "SYMBOL") {
  count_cols <- setdiff(colnames(df)[sapply(df, is.numeric)], c("GeneID", "row_sum"))
  df$median_expr <- apply(df[, count_cols, drop = FALSE], 1, median, na.rm = TRUE)
  df <- df %>%
    group_by(.data[[sym_col]]) %>%
    arrange(desc(median_expr), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    select(-median_expr)
  return(df)
}

# ==================================================================
# SECTION 1: GSE94744 (Entrez → Symbol)
# ==================================================================
file <- "GSE94744_raw_counts_GRCh38.p13_NCBI.tsv.gz"
df <- fread(file, data.table = FALSE)

cat("Original GSE94744 dimension:", dim(df), "\n")

entrez_ids <- as.character(df$GeneID)

# Report multi-mappings
all_mappings <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "SYMBOL",
                       keytype = "ENTREZID", multiVals = "list")
multi_map_ids <- names(all_mappings)[lengths(all_mappings) > 1]
cat("Entrez IDs with multiple SYMBOL mappings:", length(multi_map_ids), "\n")
if (length(multi_map_ids) > 0) {
  cat("First 5 examples:\n"); print(head(all_mappings[multi_map_ids], 5))
}

# Reproducible mapping (lex smallest symbol)
gene_symbols <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "SYMBOL",
                       keytype = "ENTREZID", multiVals = function(x) sort(x)[1])

df$SYMBOL <- gene_symbols
df <- df[!is.na(df$SYMBOL), ]

# === PINNED Ensembl v112 ===
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = 112)
cat("Connected to Ensembl v112:", mart@host, "\n")

gene_info <- getBM(attributes = c("external_gene_name", "gene_biotype"),
                   filters = "external_gene_name", values = df$SYMBOL, mart = mart)

pseudo_symbols <- gene_info$external_gene_name[
  grepl("pseudogene", gene_info$gene_biotype, ignore.case = TRUE)
]
df <- df[!(df$SYMBOL %in% pseudo_symbols), ]

# Clean zero-expression genes
count_cols <- setdiff(colnames(df)[sapply(df, is.numeric)], "GeneID")
df$row_sum <- rowSums(df[, count_cols, drop = FALSE], na.rm = TRUE)
df <- df[df$row_sum > 0, ]

# Duplicate handling with median (fixed)
df_clean <- dedup_by_median(df, "SYMBOL")

df_clean <- df_clean[, c("SYMBOL", setdiff(colnames(df_clean), c("SYMBOL", "GeneID", "row_sum")))]
write.csv(df_clean, "GSE94744_SYMBOL_clean.csv", row.names = FALSE, quote = TRUE)
cat("GSE94744 final dimension:", dim(df_clean), "\n")

# ==================================================================
# SECTION 2: counts_matrix.csv (already gene symbols)
# ==================================================================
file <- "counts_matrix.csv"
df <- read.csv(file, stringsAsFactors = FALSE)
cat("Original counts_matrix dimension:", dim(df), "\n")

gene_symbols <- as.character(df[, 1])

# PINNED Ensembl v112
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = 112)
gene_info <- getBM(attributes = c("external_gene_name", "gene_biotype"),
                   filters = "external_gene_name", values = gene_symbols, mart = mart)

pseudo_symbols <- gene_info$external_gene_name[
  grepl("pseudogene", gene_info$gene_biotype, ignore.case = TRUE)
]
df <- df[!(df[, 1] %in% pseudo_symbols), ]
df <- df[!is.na(df[, 1]) & df[, 1] != "", ]

df_clean <- dedup_by_median(df, colnames(df)[1])

colnames(df_clean)[1] <- "SYMBOL"
write.csv(df_clean, "counts_matrix_clean_no_pseudogenes.csv", row.names = FALSE, quote = TRUE)
cat("counts_matrix final dimension:", dim(df_clean), "\n")

# ==================================================================
# SECTION 3: GSE252323 (TXT.GZ)
# ==================================================================
file <- "GSE252323_LEAT_RNAseq_rawcount.txt.gz"
df <- fread(file, data.table = FALSE)
colnames(df)[1] <- "SYMBOL"
cat("Original GSE252323 dimension:", dim(df), "\n")

# PINNED Ensembl v112
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = 112)
gene_info <- getBM(attributes = c("external_gene_name", "gene_biotype"),
                   filters = "external_gene_name", values = df$SYMBOL, mart = mart)

df <- df[df$SYMBOL %in% gene_info$external_gene_name, ]
pseudo_symbols <- gene_info$external_gene_name[
  grepl("pseudogene", gene_info$gene_biotype, ignore.case = TRUE)
]
df <- df[!(df$SYMBOL %in% pseudo_symbols), ]
df <- df[!is.na(df$SYMBOL) & df$SYMBOL != "", ]

df_clean <- dedup_by_median(df, "SYMBOL")

write.csv(df_clean, "GSE252323_LEAT_RNAseq_SYMBOL_clean.csv", row.names = FALSE, quote = TRUE)
cat("GSE252323 final dimension:", dim(df_clean), "\n")

# ==================================================================
# Duplication check across all cleaned files
# ==================================================================
ds_files <- c("GSE94744_SYMBOL_clean.csv",
              "GSE134697_SYMBOL_clean.csv",
              "GSE186334_SYMBOL_clean.csv",
              "counts_matrix_clean_no_pseudogenes.csv",
              "GSE213488_SYMBOL_clean.csv",
              "GSE256068_SYMBOL_clean.csv",
              "GSE252323_LEAT_RNAseq_SYMBOL_clean.csv")

for (f in ds_files) {
  ds <- read.csv(f, stringsAsFactors = FALSE)
  cat("Duplicates in", f, ":", sum(duplicated(ds$SYMBOL)), "\n")
}

# ==================================================================
# MERGING (Intersection kept exactly as you wanted)
# ==================================================================
files <- c(
  "GSE186334_SYMBOL_clean.csv",
  "GSE134697_SYMBOL_clean.csv",
  "counts_matrix_clean_no_pseudogenes.csv",
  "GSE213488_SYMBOL_clean.csv",
  "GSE256068_SYMBOL_clean.csv",
  "GSE252323_LEAT_RNAseq_SYMBOL_clean.csv",
  "GSE94744_SYMBOL_clean.csv"
)

# Only read files that actually exist
existing_files <- files[file.exists(files)]
if (length(existing_files) < 2) stop("Need at least 2 datasets to merge. Missing files: ", 
                                     paste(files[!file.exists(files)], collapse=", "))
cat("Loading", length(existing_files), "existing files:\n", paste(existing_files, collapse="\n"), "\n")

datasets <- lapply(existing_files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  colnames(df)[1] <- "SYMBOL"
  df
})

common_genes <- Reduce(intersect, lapply(datasets, function(df) df$SYMBOL))
cat("Common genes across all 7 datasets:", length(common_genes), "\n")

datasets_common <- lapply(datasets, function(df) df[df$SYMBOL %in% common_genes, ])

merged <- Reduce(function(x, y) merge(x, y, by = "SYMBOL", all = FALSE), datasets_common)

write.csv(merged, "merged_common_genes_matrix.csv", row.names = FALSE, quote = TRUE)
cat("Merged dimension (intersection):", dim(merged), "\n")

# ==================================================================
# QC + FILTERING (Hard floor only — scientifically defensible)
# ==================================================================
merged_matrix <- read.csv("merged_common_genes_matrix.csv", check.names = FALSE)
gene_symbols <- merged_matrix[, 1]
merged_matrix <- merged_matrix[, -1]
merged_matrix <- as.data.frame(lapply(merged_matrix, as.numeric))
rownames(merged_matrix) <- gene_symbols

# QC Metrics
qc_metrics <- data.frame(
  Sample = colnames(merged_matrix),
  LibrarySize = colSums(merged_matrix, na.rm = TRUE),
  DetectedGenes = colSums(merged_matrix > 0, na.rm = TRUE),
  ZeroFraction = colSums(merged_matrix == 0, na.rm = TRUE) / nrow(merged_matrix),
  stringsAsFactors = FALSE
)

# ==================================================================
# Hard floor QC — clean, simple, fully defensible

merged_matrix <- read.csv("merged_common_genes_matrix.csv", check.names = FALSE)
gene_symbols <- merged_matrix[, 1]
merged_matrix <- merged_matrix[, -1]
merged_matrix <- as.data.frame(lapply(merged_matrix, as.numeric))
rownames(merged_matrix) <- gene_symbols

# QC Metrics
# QC Metrics
qc_metrics <- data.frame(
  Sample = colnames(merged_matrix),
  LibrarySize = colSums(merged_matrix, na.rm = TRUE),
  DetectedGenes = colSums(merged_matrix > 0, na.rm = TRUE),
  ZeroFraction = colSums(merged_matrix == 0, na.rm = TRUE) / nrow(merged_matrix),
  stringsAsFactors = FALSE
)

# ==================================================================
# Hard floor QC — clean, simple, fully defensible
# ==================================================================
qc_metrics$QC_Fail <- (qc_metrics$LibrarySize < 500000) |
  (qc_metrics$ZeroFraction > 0.75)

remove_ids <- qc_metrics$Sample[qc_metrics$QC_Fail]
cat("Samples removed by QC:", length(remove_ids), "\n")
cat("Removed sample details:\n")
print(qc_metrics[qc_metrics$QC_Fail == TRUE, ])

merged_matrix_filtered <- merged_matrix[, !colnames(merged_matrix) %in% remove_ids]

# ==================================================================
# Gold-standard low-expression filtering
# ==================================================================
# Use merged_matrix (the object that actually exists)
cpm_mat <- cpm(merged_matrix_filtered, log = FALSE)
min_samples <- ceiling(0.10 * ncol(cpm_mat))   # at least 10% of samples

keep_genes <- rowSums(cpm_mat >= 1) >= min_samples

merged_matrix_filtered <- merged_matrix_filtered[keep_genes, ]
cat("Genes retained after CPM filter:", sum(keep_genes), "\n")
# Final output
df_out <- data.frame(SYMBOL = rownames(merged_matrix_filtered),
                     merged_matrix_filtered, check.names = FALSE)
dim(df_out)
write.csv(df_out, "final_clean_matrix.csv", row.names = FALSE, quote = TRUE)
# ==================================================================
# VERIFICATION
# ==================================================================
cat("Final clean matrix dimension:", dim(df_out), "\n")
cat("Original samples:", ncol(merged_matrix), "\n")
cat("Retained samples:", ncol(merged_matrix_filtered), "\n")

reloaded <- read.csv("final_clean_matrix.csv", check.names = FALSE)
cat("Total expression matrix samples (correct):", ncol(reloaded) - 1, "\n")













