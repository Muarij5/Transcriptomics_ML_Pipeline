if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("GEOquery", quietly = TRUE))
  BiocManager::install("GEOquery")
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")
if (!requireNamespace("sva", quietly = TRUE))
  BiocManager::install("sva")
install.packages("ggVennDiagram", dependencies = TRUE)
remove.packages("ggVennDiagram")
install.packages("ggVennDiagram", dependencies = TRUE)
setwd("D:/Epilepsy_merged")

# Step 3: Load it
library(ggVennDiagram)
library(GEOquery)
library(DESeq2)
library(sva) # for ComBat
library(ggplot2)
library(dplyr)
library(lme4)
library(ggVennDiagram)
library(httr)
library(readr)
library(dplyr)
library(org.Hs.eg.db)
# ============================================
# LOAD EXPRESSION MATRIX
# ============================================
expr <- read.csv("final_clean_matrix.csv", check.names = FALSE)
gene_symbols <- expr[, 1]
expr <- expr[, -1]
rownames(expr) <- gene_symbols
expr_samples <- colnames(expr)   # vector of sample IDs
cat("Total samples in expression matrix:", length(expr_samples), "\n")
#_____________________________________
# Dataset 01
#_____________________________________
options(timeout = 3600)
gse1 <- getGEO("GSE186334",GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE)
gse1 <- gse1[[1]]
meta_raw1 <- pData(gse1)
print(colnames(meta_raw1))
# ============================================
# SELECT, RENAME AND CLEAN
# ============================================
meta_clean <- meta_raw1 %>%
  select(
    sample_id          = geo_accession,
    tissue             = source_name_ch1,
    condition          = description,
    platform           = platform_id
  ) %>%
  mutate(
    # --- Tissue ---
    tissue = recode(tissue,
                    "Cortex"      = "Cortex",
                    "Hippocampus" = "Hippocampus"
    ),
    # --- Condition / Disease group ---
    condition = recode(condition,
                       "Control"     = "Control",
                       "mTLE non-HS" = "mTLE_nonHS",
                       "mTLE+HS"     = "mTLE_HS"
    ),
    # --- Platform ---
    platform = recode(platform,
                      "GPL20301" = "Illumina_HiSeq4000"
    ),
    # --- Dataset tag ---
    dataset_id = "GSE186334",
  )
# ============================================
# VERIFY
# ============================================
cat("=== FINAL COLUMN NAMES ===\n")
print(colnames(meta_clean))

cat("\n=== DIMENSIONS ===\n")
print(dim(meta_clean))

cat("\n=== CONDITION TABLE ===\n")
print(table(meta_clean$condition))

cat("\n=== TISSUE TABLE ===\n")
print(table(meta_clean$tissue))

cat("\n=== MISSING VALUES ===\n")
print(colSums(is.na(meta_clean)))

cat("\n=== PREVIEW ===\n")
print(head(meta_clean))

# ============================================
# SAVE
# ============================================
write.csv(meta_clean, "metadata_GSE186334.csv", row.names = FALSE)
cat("\nSaved: metadata_GSE186334.csv\n")
# ============================================
# SAVE
# ============================================
write.csv(meta_clean, "metadata_GSE186334.csv", row.names = FALSE)
cat("\nSaved: metadata_GSE186334.csv\n")
#_____________________________________________________________________________________
#__DATASETS 02_______________________________________________________________________
options(timeout = 3600)
gse2 <- getGEO("GSE94744",GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE)
gse2 <- gse2[[1]]
meta_raw <- pData(gse2)
print(colnames(meta_raw))
cat("\n-- title --\n");                print(unique(meta_raw$title))
cat("\n-- source_name_ch1 --\n");     print(unique(meta_raw$source_name_ch1))
cat("\n-- characteristics_ch1 --\n"); print(unique(meta_raw$characteristics_ch1))
cat("\n-- characteristics_ch1.1 --\n"); print(unique(meta_raw$characteristics_ch1.1))
cat("\n-- molecule_ch1 --\n");        print(unique(meta_raw$molecule_ch1))
cat("\n-- instrument_model --\n");    print(unique(meta_raw$instrument_model))
cat("\n-- platform_id --\n");         print(unique(meta_raw$platform_id))
cat("\n-- library_strategy --\n");    print(unique(meta_raw$library_strategy))
cat("\n-- patient:ch1 --\n");         print(unique(meta_raw$`patient:ch1`))
cat("\n-- tissue:ch1 --\n");          print(unique(meta_raw$`tissue:ch1`))

meta_clean2 <- meta_raw %>%
  select(
    sample_id          = geo_accession,
    tissue             = `tissue:ch1`,
    platform           = platform_id
  ) %>%
  mutate(
    # --- Tissue (hippocampal subregions) ---
    tissue = recode(tissue,
                    "CA1"           = "CA1",
                    "CA3"           = "CA3",
                    "Dentate Gyrus" = "Dentate_Gyrus",
                    "Subiculum"     = "Subiculum"
    ),
    
    # --- All samples are mTLE patients, no controls ---
    condition = "mTLE",
    
    # --- Platform ---
    platform = recode(platform,
                      "GPL18460" = "Illumina_HiSeq1500"
    ),
    
    # --- Dataset tag ---
    dataset_id       = "GSE94744"
  )
# ============================================
# VERIFY
# ============================================
cat("=== DIMENSIONS ===\n");         print(dim(meta_clean2))
cat("\n=== TISSUE TABLE ===\n");     print(table(meta_clean2$tissue))
cat("\n=== CONDITION TABLE ===\n");  print(table(meta_clean2$condition))
cat("\n=== PATIENT TABLE ===\n");    print(table(meta_clean2$patient_id))
cat("\n=== MISSING VALUES ===\n");   print(colSums(is.na(meta_clean2)))
cat("\n=== PREVIEW ===\n");          print(head(meta_clean2))

# ============================================
# SAVE
# ============================================
write.csv(meta_clean2, "metadata_GSE94744.csv", row.names = FALSE)
cat("\nSaved: metadata_GSE94744.csv\n")
#___________________________________
#DATASET 03
#_____________________________________
options(timeout = 3600)
gse3 <- getGEO("GSE134697", GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE)
gse3 <- gse3[[1]]
meta_raw3 <- pData(gse3)
print(colnames(meta_raw3))
cat("\n-- title --\n");                  print(unique(meta_raw3$title))
cat("\n-- source_name_ch1 --\n");        print(unique(meta_raw3$source_name_ch1))
cat("\n-- description --\n");            print(unique(meta_raw3$description))
cat("\n-- characteristics_ch1 --\n");    print(unique(meta_raw3$characteristics_ch1))
cat("\n-- characteristics_ch1.1 --\n");  print(unique(meta_raw3$characteristics_ch1.1))
cat("\n-- characteristics_ch1.2 --\n");  print(unique(meta_raw3$characteristics_ch1.2))
cat("\n-- molecule_ch1 --\n");           print(unique(meta_raw3$molecule_ch1))
cat("\n-- instrument_model --\n");       print(unique(meta_raw3$instrument_model))
cat("\n-- platform_id --\n");            print(unique(meta_raw3$platform_id))
cat("\n-- library_strategy --\n");       print(unique(meta_raw3$library_strategy))
cat("\n-- condition:ch1 --\n");          print(unique(meta_raw3$`condition:ch1`))
cat("\n-- individual:ch1 --\n");         print(unique(meta_raw3$`individual:ch1`))
cat("\n-- tissue:ch1 --\n");             print(unique(meta_raw3$`tissue:ch1`))

meta_clean3 <- meta_raw3 %>%
  select(
    sample_id           = geo_accession,
    tissue              = `tissue:ch1`,
    condition           = `condition:ch1`,
    platform            = platform_id
  ) %>%
  mutate(
    # --- Tissue ---
    tissue = recode(tissue,
                    "Neocortex"   = "Neocortex",
                    "Hippocampus" = "Hippocampus"
    ),
    
    # --- Condition ---
    condition = recode(condition,
                       "mesial temporal lobe epilepsy" = "mTLE",
                       "without epilepsy"              = "Control"
    ),
    
    # --- Platform ---
    platform = recode(platform,
                      "GPL16791" = "Illumina_HiSeq2500"
    ),
    # --- Dataset tag ---
    dataset_id       = "GSE134697"
  )
# ============================================
# VERIFY
# ============================================
cat("=== DIMENSIONS ===\n");        print(dim(meta_clean3))
cat("\n=== TISSUE TABLE ===\n");    print(table(meta_clean3$tissue))
cat("\n=== CONDITION TABLE ===\n"); print(table(meta_clean3$condition))
cat("\n=== MISSING VALUES ===\n");  print(colSums(is.na(meta_clean3)))
cat("\n=== PREVIEW ===\n");         print(head(meta_clean3))
# ============================================
# SAVE
# ============================================
write.csv(meta_clean3, "metadata_GSE134697.csv", row.names = FALSE)
cat("\nSaved: metadata_GSE134697.csv\n")
#__________________________
# Dataset 04
#__________________________
options(timeout = 3600)
gse4 <- getGEO("GSE213488", AnnotGPL = FALSE, 
               getGPL = FALSE)
gse4 <- gse4[[1]]
meta_raw4 <- pData(gse4)
print(colnames(meta_raw4))

# ============================================================
# INSPECT ALL RELEVANT COLUMNS — GSE213488
# ============================================================

cat("\n-- title --\n");                  print(unique(meta_raw4$title))
cat("\n-- source_name_ch1 --\n");        print(unique(meta_raw4$source_name_ch1))
cat("\n-- description --\n");            print(unique(meta_raw4$description))
cat("\n-- characteristics_ch1 --\n");    print(unique(meta_raw4$characteristics_ch1))
cat("\n-- characteristics_ch1.1 --\n");  print(unique(meta_raw4$characteristics_ch1.1))
cat("\n-- molecule_ch1 --\n");           print(unique(meta_raw4$molecule_ch1))
cat("\n-- instrument_model --\n");       print(unique(meta_raw4$instrument_model))
cat("\n-- platform_id --\n");            print(unique(meta_raw4$platform_id))
cat("\n-- library_strategy --\n");       print(unique(meta_raw4$library_strategy))
cat("\n-- disease state:ch1 --\n");      print(unique(meta_raw4$`disease state:ch1`))
cat("\n-- tissue:ch1 --\n");             print(unique(meta_raw4$`tissue:ch1`))

meta_clean4 <- meta_raw4 %>%
  select(
    sample_id           = geo_accession,
    tissue              = source_name_ch1,
    condition           = `disease state:ch1`,
    platform            = platform_id,
    title               = title
  ) %>%
  mutate(
    
    # --- Tissue ---
    tissue = recode(tissue,
                    "Frontal Lobe" = "Frontal_Lobe"
    ),
    
    # --- Condition ---
    condition = recode(condition,
                       "FCD IIa"           = "FCD_IIa",
                       "FCD IIb"           = "FCD_IIb",
                       "Control (autopsy)" = "Control"
    ),
    
    # --- Platform ---
    platform            = recode(platform, "GPL16791" = "Illumina_HiSeq2500"),
    # --- Dataset tag ---
    dataset_id       = "GSE213488"
  ) %>%
  select(-title)  # drop title after extracting matter

# ============================================
# VERIFY
# ============================================
cat("=== DIMENSIONS ===\n");         print(dim(meta_clean4))
cat("\n=== CONDITION TABLE ===\n");  print(table(meta_clean4$condition))
cat("\n=== MATTER TABLE ===\n");     print(table(meta_clean4$matter))
cat("\n=== MISSING VALUES ===\n");   print(colSums(is.na(meta_clean4)))
cat("\n=== PREVIEW ===\n");          print(head(meta_clean4))

# ============================================
# SAVE
# ============================================
write.csv(meta_clean4, "metadata_GSE213488.csv", row.names = FALSE)
cat("\nSaved: metadata_GSE213488.csv\n")
#_____________________________________
#Dataset 05
#______________________________________
options(timeout = 3600)
gse5 <- getGEO("GSE256068", GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE)
gse5 <- gse5[[1]]
meta_raw5 <- pData(gse5)
print(colnames(meta_raw5))
cat("\n-- title --\n");                  print(unique(meta_raw5$title))
cat("\n-- source_name_ch1 --\n");        print(unique(meta_raw5$source_name_ch1))
cat("\n-- description --\n");            print(unique(meta_raw5$description))
cat("\n-- characteristics_ch1 --\n");    print(unique(meta_raw5$characteristics_ch1))
cat("\n-- characteristics_ch1.1 --\n");  print(unique(meta_raw5$characteristics_ch1.1))
cat("\n-- characteristics_ch1.2 --\n");  print(unique(meta_raw5$characteristics_ch1.2))
cat("\n-- characteristics_ch1.3 --\n");  print(unique(meta_raw5$characteristics_ch1.3))
cat("\n-- characteristics_ch1.4 --\n");  print(unique(meta_raw5$characteristics_ch1.4))
cat("\n-- molecule_ch1 --\n");           print(unique(meta_raw5$molecule_ch1))
cat("\n-- instrument_model --\n");       print(unique(meta_raw5$instrument_model))
cat("\n-- platform_id --\n");            print(unique(meta_raw5$platform_id))
cat("\n-- library_strategy --\n");       print(unique(meta_raw5$library_strategy))
cat("\n-- age:ch1 --\n");               print(unique(meta_raw5$`age:ch1`))
cat("\n-- Sex:ch1 --\n");               print(unique(meta_raw5$`Sex:ch1`))
cat("\n-- disease state:ch1 --\n");     print(unique(meta_raw5$`disease state:ch1`))
cat("\n-- seisure frequency:ch1 --\n"); print(unique(meta_raw5$`seisure frequency:ch1`))
cat("\n-- tissue:ch1 --\n");            print(unique(meta_raw5$`tissue:ch1`))
library(dplyr)

meta_clean5 <- meta_raw5 %>%
  select(
    sample_id           = geo_accession,
    tissue              = `tissue:ch1`,
    Gender                 = `Sex:ch1`,
    age                 = `age:ch1`,
    condition           = `disease state:ch1`,
    platform            = platform_id
  ) %>%
  mutate(
    # --- Tissue ---
    tissue = recode(tissue,
                    "Temporal"    = "Temporal",
                    "Frontal"     = "Frontal",
                    "Parietal"    = "Parietal",
                    "Hippocampus" = "Hippocampus",
                    "Occipital"   = "Occipital"
    ),
    # --- Condition ---
    condition = recode(condition,
                       "TSC"                = "TSC",
                       "FCD2a"              = "FCD_IIa",
                       "FCD2b"              = "FCD_IIb",
                       "TLE-HS"             = "mTLE_HS",
                       "ControlCortex"      = "Control",
                       "ControlHippocampus" = "Control"
    ),
    
    # --- Platform ---
    platform            = recode(platform, "GPL24676" = "Illumina_NovaSeq6000"),
    
    
    # --- Dataset tag ---
    dataset_id       = "GSE256068"
  )

# ============================================
# VERIFY
# ============================================
cat("=== DIMENSIONS ===\n");               print(dim(meta_clean5))
cat("\n=== CONDITION TABLE ===\n");         print(table(meta_clean5$condition))
cat("\n=== TISSUE TABLE ===\n");            print(table(meta_clean5$tissue))
cat("\n=== SEX TABLE ===\n");               print(table(meta_clean5$sex))
cat("\n=== AGE SUMMARY ===\n");             print(summary(meta_clean5$age))
cat("\n=== SEIZURE FREQ SUMMARY ===\n");    print(summary(meta_clean5$seizure_frequency))
cat("\n=== MISSING VALUES ===\n");          print(colSums(is.na(meta_clean5)))
cat("\n=== PREVIEW ===\n");                 print(head(meta_clean5))

# ============================================
# SAVE
# ============================================
write.csv(meta_clean5, "metadata_GSE256068.csv", row.names = FALSE)
cat("\nSaved: metadata_GSE256068.csv\n")
#_____________________________________
# Dataset 06
#_____________________________________
expr_matrix <- read.csv("final_clean_matrix.csv", row.names = 1)
expr_samples <- colnames(expr_matrix)

options(timeout = 3600)
gse6 <- getGEO("GSE252323", GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE)
gse6 <- gse6[[1]]
meta_raw6 <- pData(gse6)

print(colnames(meta_raw6))
cat("\n-- title --\n");               print(unique(meta_raw6$title))
cat("\n-- source_name_ch1 --\n");     print(unique(meta_raw6$source_name_ch1))
cat("\n-- description --\n");         print(unique(meta_raw6$description))
cat("\n-- characteristics_ch1 --\n"); print(unique(meta_raw6$characteristics_ch1))
cat("\n-- molecule_ch1 --\n");        print(unique(meta_raw6$molecule_ch1))
cat("\n-- instrument_model --\n");    print(unique(meta_raw6$instrument_model))
cat("\n-- platform_id --\n");         print(unique(meta_raw6$platform_id))
cat("\n-- library_strategy --\n");    print(unique(meta_raw6$library_strategy))
cat("\n-- tissue:ch1 --\n");          print(unique(meta_raw6$`tissue:ch1`))

# ============================================
# STEP 1: CHECK TITLE FORMATS
# ============================================
gse252323_titles <- c("Control_C01", "Control_C02", "Control_C03", "Control_C04",
                      "LEAT_P01", "LEAT_P02", "LEAT_P03", "LEAT_P04", "LEAT_P05",
                      "LEAT_P06", "LEAT_P07", "LEAT_P08", "LEAT_P09", "LEAT_P10",
                      "LEAT_P11", "LEAT_P12", "LEAT_P13", "LEAT_P14", "LEAT_P15",
                      "LEAT_P16", "LEAT_P17", "LEAT_P18", "LEAT_P19", "LEAT_P20",
                      "LEAT_P21")

cat("Matching with full titles (Control_C01 format):\n")
print(sum(gse252323_titles %in% expr_samples))

plain_titles <- gsub("Control_|LEAT_", "", gse252323_titles)
cat("\nMatching with plain titles (C01, P01 format):\n")
print(sum(plain_titles %in% expr_samples))

# ============================================
# STEP 2: REBUILD METADATA - KEEP ONLY CONTROLS
# ============================================
meta_clean6 <- meta_raw6 %>%
  mutate(
    sample_id = sub("^[^_]+_", "", title),  # "Control_C01" → "C01"
    condition = case_when(
      grepl("Control", title, ignore.case = TRUE) ~ "Control",
      grepl("LEAT",    title, ignore.case = TRUE) ~ "BRAF_V600E_LEAT",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(condition == "Control") %>%        # Keep ONLY Controls
  select(
    sample_id,
    geo_accession,
    condition,
    platform = platform_id
  ) %>%
  mutate(
    tissue     = "Surgical_Brain",
    platform   = "Illumina_NovaSeq6000",
    dataset_id = "GSE252323"
  )

# VERIFY
cat("Control samples found:", nrow(meta_clean6), "\n")
cat("sample_ids extracted:", paste(meta_clean6$sample_id, collapse = ", "), "\n")
cat("Matched:", sum(meta_clean6$sample_id %in% expr_samples), "/", nrow(meta_clean6), "\n")

write.csv(meta_clean6, "metadata_GSE252323.csv", row.names = FALSE)
cat("Saved: metadata_GSE252323.csv\n")

#Dataset 07
#______________________________
options(timeout = 3600)
gse7 <- getGEO("GSE310367", GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE)
gse7 <- gse7[[1]]
meta_raw7 <- pData(gse7)

print(colnames(meta_raw7))
cat("\n-- title --\n");                  print(unique(meta_raw7$title))
cat("\n-- source_name_ch1 --\n");        print(unique(meta_raw7$source_name_ch1))
cat("\n-- description --\n");            print(unique(meta_raw7$description))
cat("\n-- characteristics_ch1 --\n");    print(unique(meta_raw7$characteristics_ch1))
cat("\n-- characteristics_ch1.1 --\n");  print(unique(meta_raw7$characteristics_ch1.1))
cat("\n-- molecule_ch1 --\n");           print(unique(meta_raw7$molecule_ch1))
cat("\n-- instrument_model --\n");       print(unique(meta_raw7$instrument_model))
cat("\n-- platform_id --\n");            print(unique(meta_raw7$platform_id))
cat("\n-- library_strategy --\n");       print(unique(meta_raw7$library_strategy))
cat("\n-- group:ch1 --\n");              print(unique(meta_raw7$`group:ch1`))
cat("\n-- tissue:ch1 --\n");             print(unique(meta_raw7$`tissue:ch1`))

library(dplyr)

meta_clean7 <- meta_raw7 %>%
  select(
    sample_id           = geo_accession,
    tissue              = `tissue:ch1`,
    condition           = `group:ch1`,
    platform            = platform_id
  ) %>%
  mutate(
    # --- Tissue ---
    tissue = recode(tissue,
                    "brain" = "Brain"
    ),
    
    # --- Condition ---
    condition = recode(condition,
                       "Control" = "Control",
                       "PWE"     = "PWE",
                       "RE"      = "Rasmussen_Encephalitis"
    ),
    
    # --- Platform ---
    platform            = recode(platform, "GPL24676" = "Illumina_NovaSeq6000"),
    # --- Dataset tag ---
    dataset_id       = "GSE310367"
  )

# ============================================
# VERIFY
# ============================================
cat("=== DIMENSIONS ===\n");        print(dim(meta_clean7))
cat("\n=== CONDITION TABLE ===\n"); print(table(meta_clean7$condition))
cat("\n=== MISSING VALUES ===\n");  print(colSums(is.na(meta_clean7)))
cat("\n=== PREVIEW ===\n");         print(head(meta_clean7))

# ============================================
# SAVE
# ============================================
write.csv(meta_clean7, "metadata_GSE310367.csv", row.names = FALSE)
cat("\nSaved: metadata_GSE310367.csv\n")

# LOAD ALL 7 METADATA FILES
# ============================================
meta1 <- read.csv("metadata_GSE186334.csv")
meta2 <- read.csv("metadata_GSE94744.csv")
meta3 <- read.csv("metadata_GSE134697.csv")
meta4 <- read.csv("metadata_GSE213488.csv")
meta5 <- read.csv("metadata_GSE252323.csv")
meta6 <- read.csv("metadata_GSE256068.csv")
meta7 <- read.csv("metadata_GSE310367.csv")

expr_samples <- read.csv("final_clean_matrix.csv")

meta_list <- list(
  GSE186334 = meta1,
  GSE94744  = meta2,
  GSE134697 = meta3,
  GSE213488 = meta4,
  GSE252323 = meta5,
  GSE256068 = meta6,
  GSE310367 = meta7
)

for (name in names(meta_list)) {
  meta      <- meta_list[[name]]
  meta_ids  <- meta$sample_id
  
  matched   <- sum(meta_ids %in% expr_samples)
  unmatched <- sum(!meta_ids %in% expr_samples)
  
  cat("===================\n")
  cat("Dataset:", name, "\n")
  cat("Samples in metadata:         ", length(meta_ids), "\n")
  cat("Matched to expression matrix:", matched, "\n")
  cat("In metadata NOT in matrix:   ", unmatched, "\n")
  
  if (unmatched > 0) {
    cat("Unmatched sample IDs:\n")
    print(meta_ids[!meta_ids %in% expr_samples])
  }
}

# ============================================
# MERGE & SAVE METADATA
# ============================================
final_cols <- c("sample_id", "tissue", "condition", "platform", "dataset_id")

merged_meta <- do.call(rbind, lapply(meta_list, function(m) {
  m[, intersect(final_cols, colnames(m)), drop = FALSE]
}))

write.csv(merged_meta, "merged_metadata.csv", row.names = FALSE)
cat("Merged metadata saved to merged_metadata.csv\n")
cat("Dimensions:", nrow(merged_meta), "rows x", ncol(merged_meta), "cols\n")
merged_meta=read.csv("merged_metadata.csv")
dim(merged_meta)
# ============================================
# OVERALL SUMMARY
# ============================================
all_meta_ids <- unlist(lapply(meta_list, function(m) m$sample_id))
cat("\n===================\n")
cat("Total metadata samples across all datasets:", length(all_meta_ids), "\n")
cat("Total expression matrix samples:           ", length(expr_samples), "\n")
cat("Overall matched:                           ", sum(all_meta_ids %in% expr_samples), "\n")

cat("First 10 expression matrix sample names:\n")
print(head(expr_samples, 10))

cat("\nFirst 10 GSE252323 metadata sample IDs:\n")
print(head(meta5$sample_id, 10))

cat("\n=== SAMPLE NAME FORMAT CHECK ===\n")
for (name in names(meta_list)) {
  cat("\n", name, "- first 3 metadata IDs:\n")
  print(head(meta_list[[name]]$sample_id, 3))
}