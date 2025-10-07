###############################################################
# Script: Exploratory Data Analysis (EDA)
# Author: Prasad Chaskar
# Date: 07-10-25
# Description:
# - Perform PCA and hierarchical clustering heatmap
# - Explore relationships and sample grouping in gene expression data
###############################################################

# ---------------------------
# 0. Setup & environment
# ---------------------------

library(future)
options(future.globals.maxSize = 100000 * 1024^2)  # Allow large global objects
options(oligoOptions = list(nCores = 12))           # Control parallel threads

# plan("multisession", workers = 32)               # Optional: parallel computation

# ---------------------------
# 1. Working directory setup
# ---------------------------

library(rstudioapi)
script_path <- getActiveDocumentContext()$path      # Locate current script path
dirpath <- dirname(script_path)
setwd(dirname(dirpath))                             # Move one directory up
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
  library(pd.huex.1.0.st.v2)                        # Platform definition
  library(huex10sttranscriptcluster.db)             # Annotation database
  library(AnnotationDbi)
  library(nclust)                                   # For coldmap heatmaps
  library(PCAtools)                                 # For PCA visualization
})

# ---------------------------
# 3. Load normalized expression data
# ---------------------------

final_eset <- readRDS("./RDS/GSE224644_RMA_expression_gene_symbols.rds")
cat(
  "ExpressionSet loaded with",
  nrow(final_eset),
  "features and",
  ncol(final_eset),
  "samples.\n"
)

# ---------------------------
# 4. Extract normalized expression data
# ---------------------------

ex <- exprs(final_eset)     # log2 RMA-normalized expression values
cat("Expression matrix dimensions (features × samples):",
    dim(ex),
    "\n")
head(ex[, 1:5])             # Preview the first few samples

# ---------------------------
# 5. Extract phenotype metadata
# ---------------------------

metadata <- pData(final_eset)     # GEO sample metadata
cat("Metadata dimensions:", dim(metadata), "\n")
metadata[1:5, 1:5]

# ---------------------------
# 1. Parse 'title' into separate fields
# ---------------------------

# Example: "prostate_cancer_primary_1 (STAMPEDE)"
metadata <- metadata %>%
  mutate(
    sample_label = sub(" \\(.*\\)$", "", title),
    # remove "(STAMPEDE)"
    study_id = sub(".*\\((.*)\\)$", "\\1", title)           # extract content inside parentheses
  )

# ---------------------------
# 2. Parse 'characteristics_ch1' into multiple columns
# ---------------------------

# Example: "gleason score: 9; age at randomization: 66; psa at randomization (ng/ml): 31; tissue: FFPE prostate cancer"

metadata <- metadata %>%
  mutate(characteristics_ch1 = as.character(characteristics_ch1)) %>%
  separate_wider_delim(
    characteristics_ch1,
    delim = ";",
    names = c("gleason", "age", "psa", "tissue"),
    too_few = "align_start"
  ) %>%
  mutate(
    gleason = trimws(sub("gleason score:\\s*", "", gleason)),
    age = as.numeric(trimws(sub(
      "age at randomization:\\s*", "", age
    ))),
    psa = as.numeric(trimws(
      sub("psa at randomization \\(ng/ml\\):\\s*", "", psa)
    )),
    tissue = trimws(sub("tissue:\\s*", "", tissue))
  )

# ---------------------------
# 3. Clean up and reorder
# ---------------------------
metadata <- metadata %>%
  relocate(sample_label, study_id, gleason, age, psa, tissue, .before = 1)

# ---------------------------
# 4. Check results
# ---------------------------
head(metadata[, c("sample_label", "study_id", "gleason", "age", "psa", "tissue")])

metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$geo_accession

metadata <- metadata %>%
  mutate(
    metastatic_guess = case_when(
      gleason >= 9 | psa > 100 ~ "Likely metastatic",
      gleason >= 8 | psa > 20 ~ "Possible metastatic",
      TRUE ~ "Likely non-metastatic"
    )
  )

table(metadata$metastatic_guess, useNA = "ifany")

# ---------------------------
# 6. Ensure matching samples between expression and metadata
# ---------------------------

ex <- ex[, which(colnames(ex) %in% rownames(metadata))]  # Filter only overlapping samples
stopifnot(all(colnames(ex) == rownames(metadata)))       # Validate order alignment

# ---------------------------
# 7. Hierarchical clustering and heatmap visualization
# ---------------------------

cat("Generating scaled heatmap for hierarchical clustering...\n")

# Scale expression by gene (row-wise z-score)

scaled_data <- t(scale(t(ex)))

# Perform hierarchical clustering using Ward’s method

hist <- coldmap(scaled_data, method = "ward")

saveRDS(hist, "./RDS/hist_GSE224644.rds")
saveRDS(metadata, "./RDS/metadata_GSE224644.rds")

# Visualize the clustered heatmap

coldmap(
  scaled_data,
  clust = hist,
  saturation = TRUE,
  ctag.space = 1,
  rlab = list(list(
    c("CD3[DEG]", "TP53", "PTEN", "SPOP", "AR", "ERG", "MKI67")
  )),
  # ctag = make_tag(
  # metadata,
  # varnames = c("metastatic_guess"),
  # cols = c("green4")
  # )
)

# ---------------------------
# 8. Principal Component Analysis (PCA)
# ---------------------------

cat("Performing PCA...\n")

# Remove the lowest 10% variable genes to reduce noise

p <- pca(ex, metadata = metadata, removeVar = 0.1)

# Basic PCA biplot

biplot(p)

# PCA biplot with loadings

biplot(p, showLoadings = TRUE, lab = NULL)

# Colored PCA biplot by experimental factor (e.g., agent.ch1)

p_plot <- biplot(
  p,
  colby = 'agent.ch1',
  legendPosition = 'right',
  labSize = 5,
  xlim = c(-35, 35),
  ylim = c(-35, 35),
  pointSize = 4
)

# ---------------------------
# 9. Save PCA plot to file
# ---------------------------

cat("Saving PCA plot to ./images...\n")
dir.create("./images", showWarnings = FALSE)

png(
  paste("./images/EDA_PCA", "All_Sample", Sys.Date(), "PCA.png", sep = "_"),
  width = 11.68,
  height = 8.62,
  units = 'in',
  res = 500
)
plot(p_plot)
dev.off()

cat("EDA complete. PCA and heatmap generated successfully.\n")
