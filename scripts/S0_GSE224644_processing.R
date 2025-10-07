###############################################################
# Script: RMA normalization for Human Exon 1.0 ST v2 arrays
# Author: Prasad Chaskar
# Date: 07-10-25
# Description:
#   - Decompress CEL files
#   - Perform RMA normalization in memory-safe chunks
#   - Attach phenotype data directly from GEO
#   - Annotate probes and aggregate to gene symbols
###############################################################

# ---------------------------
# 0. Setup & environment
# ---------------------------
library(future)
options(future.globals.maxSize = 100000 * 1024^2)  # Allow large global objects
options(oligoOptions = list(nCores = 12))           # Limit CPU cores to control RAM
# plan("multisession", workers = 32)               # Optional parallelization

# ---------------------------
# 1. Working directory setup
# ---------------------------
library(rstudioapi)
script_path <- getActiveDocumentContext()$path
dirpath <- dirname(script_path)
setwd(dirname(dirpath))  # Move one level up
cat("Working directory set to:", getwd(), "\n")

# ---------------------------
# 2. Load required packages
# ---------------------------
suppressPackageStartupMessages({
  library(affyio)
  library(affxparser)
  library(oligo)
  library(GEOquery)
  library(tidyverse)
  library(Biobase)
  library(pd.huex.1.0.st.v2)                       # Platform definition
  library(huex10sttranscriptcluster.db)            # Annotation database
  library(AnnotationDbi)
})

# ---------------------------
# 3. Decompress CEL.gz files (if needed)
# ---------------------------
gz_files <- list.files("../00_RAW_DATA/data/", pattern = "\\.CEL\\.gz$", full.names = TRUE)
for (f in gz_files) {
  out <- sub("\\.gz$", "", f)
  if (!file.exists(out)) {
    system2("gunzip", c("-c", shQuote(f)), stdout = out)
    cat("Decompressed:", basename(f), "\n")
  }
}

# ---------------------------
# 4. Read CEL files
# ---------------------------
celfiles <- list.celfiles("../00_RAW_DATA/data/", full.names = TRUE)
cat("Total CEL files found:", length(celfiles), "\n")

# Optional: check platform
hdr <- readCelHeader(celfiles[1])
cat("Detected chip type:", hdr$chiptype, "\n")

# ---------------------------
# 5. Chunk-wise RMA normalization
# ---------------------------
chunks <- split(celfiles, ceiling(seq_along(celfiles)/500))  # ~500 CELs per chunk
exprs_list <- list()

for (i in seq_along(chunks)) {
  cat("Processing chunk", i, "of", length(chunks), "\n")
  fs <- read.celfiles(chunks[[i]])
  eset <- rma(fs, target = "core")  # Gene-level summarization
  exprs_list[[i]] <- exprs(eset)
  rm(fs, eset); gc()
}

# Combine all chunks
exprs_mat <- do.call(cbind, exprs_list)

# Clean sample names: remove .CEL and trailing _001, _002, etc.
colnames(exprs_mat) <- gsub("\\.CEL$", "", basename(colnames(exprs_mat)))
colnames(exprs_mat) <- sub("_\\d+$", "", colnames(exprs_mat))
head(colnames(exprs_mat))

# ---------------------------
# 6. Retrieve phenotype metadata from GEO
# ---------------------------
gsm_ids <- colnames(exprs_mat)
gsm_list <- lapply(gsm_ids, getGEO)

# Flatten GSM metadata (1 row per sample)
meta_data <- lapply(gsm_list, function(gsm) {
  pdata <- Meta(gsm)
  pdata_flat <- sapply(pdata, function(x) paste(x, collapse = "; "))
  df <- as.data.frame(t(pdata_flat), stringsAsFactors = FALSE)
  rownames(df) <- gsm@header$geo_accession
  return(df)
})

# Combine into single data.frame
meta_data_combined <- bind_rows(meta_data, .id = "GSM")
rownames(meta_data_combined) <- meta_data_combined$geo_accession
meta_data_combined$GSM <- NULL

# ---------------------------
# 7. Create ExpressionSet
# ---------------------------
pheno_adf <- AnnotatedDataFrame(meta_data_combined)
final_eset <- ExpressionSet(assayData = exprs_mat, phenoData = pheno_adf)
cat("Final ExpressionSet created with", ncol(final_eset), "samples\n")

# ---------------------------
# 8. Annotate probes and aggregate to gene symbols
# ---------------------------
annot <- select(huex10sttranscriptcluster.db,
                keys = rownames(final_eset),
                columns = c("SYMBOL", "GENENAME"),
                keytype = "PROBEID")

# Keep only probes with gene symbols
annot <- annot[!is.na(annot$SYMBOL), ]
annot <- annot[match(rownames(final_eset), annot$PROBEID), ]

# Attach probe annotation to ExpressionSet
featureData(final_eset) <- AnnotatedDataFrame(annot)

# Aggregate probe-level data to gene symbols (average expression per gene)
exprs_mat <- exprs(final_eset)
gene_symbols <- annot$SYMBOL

exprs_agg <- aggregate(exprs_mat, by = list(GeneSymbol = gene_symbols), FUN = mean)
rownames(exprs_agg) <- exprs_agg$GeneSymbol
exprs_agg$GeneSymbol <- NULL

# Create new ExpressionSet with gene symbols as rownames
final_eset <- ExpressionSet(
  assayData = as.matrix(exprs_agg),
  phenoData = phenoData(final_eset)
)

# ---------------------------
# Fix annotation before using SYMBOL as rownames
# ---------------------------
gene_annot <- unique(annot[, c("SYMBOL", "GENENAME")])

# Remove entries with missing or duplicate SYMBOLs
gene_annot <- gene_annot[!is.na(gene_annot$SYMBOL) & gene_annot$SYMBOL != "", ]
gene_annot <- gene_annot[!duplicated(gene_annot$SYMBOL), ]

# Now safely assign rownames
rownames(gene_annot) <- gene_annot$SYMBOL

featureData(final_eset) <- AnnotatedDataFrame(gene_annot)

cat("Final ExpressionSet now has gene symbols as rownames.\n")

# ---------------------------
# 9. Save annotated ExpressionSet
# ---------------------------
saveRDS(final_eset, "./RDS/GSE224644_RMA_expression_gene_symbols.rds")
cat("Saved: ./RDS/GSE224644_RMA_expression_gene_symbols.rds\n")

# ---------------------------
# 10. Test: verify rownames
# ---------------------------
ex <- exprs(final_eset)
cat("Matrix dimensions:", dim(ex), "\n")
cat("Example gene symbols:\n")
print(head(rownames(ex)))
# ---------------------------
# 10. Clean up
# ---------------------------
rm(list = ls())
.rs.restartR()
