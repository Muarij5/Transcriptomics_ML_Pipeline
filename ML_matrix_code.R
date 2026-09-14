# ── REQUIRED LIBRARIES ───────────────────────────────────────
library(dplyr)
setwd("D:/Epilepsy_merged")

# Load VST matrix
ml_matrix <- read.csv("ML_matrix_VST.csv", row.names = 1, check.names = FALSE)

# Load DEG gene names
DEGs_mTLE <- rownames(read.csv("DEGs_mTLE_vs_Control_final.csv", row.names = 1))
DEGs_FCD  <- rownames(read.csv("DEGs_FCD_vs_Control_final.csv",  row.names = 1))
DEGs_TSC  <- rownames(read.csv("DEGs_TSC_vs_Control_final.csv",  row.names = 1))

# Union of all DEGs
union_genes <- unique(c(DEGs_mTLE, DEGs_FCD, DEGs_TSC))
cat("Total union DEGs:", length(union_genes), "\n")

# Subset VST matrix to only DEG columns
keep <- intersect(union_genes, colnames(ml_matrix))
cat("DEGs matched in matrix:", length(keep), "\n")
ml_final <- ml_matrix[, c(keep, "condition")]

# Rename condition to Label
colnames(ml_final)[colnames(ml_final) == "condition"] <- "Label"

# Convert Label to binary
ml_final$Label <- ifelse(ml_final$Label == "Control", 0, 1)

# Fix gene column names
colnames(ml_final) <- gsub("\\.", "-", colnames(ml_final))

# Save
write.csv(ml_final, "ML_matrix_DEGs_only.csv", row.names = TRUE)
cat("Samples:", nrow(ml_final), "\n")
cat("DEG features:", ncol(ml_final) - 1, "\n")
print(table(ml_final$Label))