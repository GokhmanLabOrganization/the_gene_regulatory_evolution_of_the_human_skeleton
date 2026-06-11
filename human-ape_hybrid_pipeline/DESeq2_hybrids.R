# This script computes p-values and ASE fold change based on DESeq2.
# See DESeq2 manual for further information
# http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html

library("DESeq2", quietly=TRUE, verbose=FALSE)
library("GenomicRanges", quietly=TRUE, verbose=FALSE)
library("qvalue", quietly=TRUE, verbose=FALSE)
library(ggplot2, quietly=TRUE, verbose=FALSE)
library(apeglm, quietly=TRUE, verbose=FALSE)

# ---------- Function to load counts ----------
load_counts <- function(species_folder, suffix) {
  files <- list.files(species_folder, pattern = paste0(suffix, "$"), full.names = TRUE)
  if (length(files) == 0) stop("No count files found in: ", species_folder)
  
  counts_list <- lapply(files, function(f) {
    dat <- read.table(f, 
                      header = TRUE, 
                      sep = "\t", 
                      stringsAsFactors = FALSE)
    # extract only gene, ref_counts, alt_counts
    dat <- dat[, c("gene", "ref_counts", "alt_counts")]
    sample_name <- sub(paste0(suffix, "$"), "", basename(f))
    # rename columns to sample-specific
    colnames(dat) <- c("gene", paste0(sample_name, "_ref"), paste0(sample_name, "_alt"))
    dat
  })
  
  # merge on gene
  counts <- Reduce(function(x, y) merge(x, y, by = "gene"), counts_list)
  
  # make gene rownames
  rownames(counts) <- counts$gene
  counts <- counts[, -1, drop = FALSE]  # drop redundant gene column (now rownames)
  
  return(counts)
}

# ---------- Function to load TPMs ----------
load_tpms <- function(species_folder, suffix) {
  files <- list.files(species_folder, pattern = paste0(suffix, "$"), full.names = TRUE)
  if (length(files) == 0) stop("No TPM-containing files found in: ", species_folder)
  
  tpm_list <- lapply(files, function(f) {
    dat <- read.table(f,
                      header = TRUE,
                      sep = "\t",
                      stringsAsFactors = FALSE)
    # Expect columns from GetGeneASEbyReads.py script:
    # gene, chrom, ref_counts, alt_counts, no_ase_counts, ambig_ase_counts,
    # ase_value, ref_TPM, alt_TPM, total_TPM
    dat <- dat[, c("gene", "ref_TPM", "alt_TPM", "total_TPM")]
    
    sample_name <- sub(paste0(suffix, "$"), "", basename(f))
    
    colnames(dat) <- c(
      "gene",
      paste0(sample_name, "_ref_tpm"),
      paste0(sample_name, "_alt_tpm"),
      paste0(sample_name, "_total_tpm")
    )
    dat
  })
  
  # merge on gene
  tpms <- Reduce(function(x, y) merge(x, y, by = "gene"), tpm_list)
  
  rownames(tpms) <- tpms$gene
  tpms <- tpms[, -1, drop = FALSE]  # drop gene column, keep as rownames
  
  return(tpms)
}


# ---------- Run DESeq2 ----------
run_deseq2 <- function(counts, coldata, test_columns) {
  # full design: all columns in coldata (except "sample")
  all_covariates <- setdiff(colnames(coldata), "sample")
  full_design <- as.formula(paste("~", paste(all_covariates, collapse = " + ")))
  
  # reduced design: all covariates except the ones to test
  reduced_covariates <- setdiff(all_covariates, test_columns)
  reduced_design <- if (length(reduced_covariates) == 0) ~1 else as.formula(paste("~", paste(reduced_covariates, collapse = " + ")))
  
  
  mm <- model.matrix(full_design, data = coldata)
  cat("Model matrix columns:", ncol(mm), "  rank:", qr(mm)$rank, "\n")
  
  # construct DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = counts, 
                                colData = coldata, 
                                design = full_design)
  
  # run DESeq with LRT
  dds <- DESeq(dds, 
               test = "LRT", 
               betaPrior=FALSE, 
               full=full_design, 
               reduced=reduced_design)
  
  print(resultsNames(dds))
  
  res <- results(dds)
  
  coef <- grep("human_vs_", resultsNames(dds)) #the coef parameter determines the column by which to shrink the data (the log-change, in my case - chi vs hum. See resultNames(dds) to see which column fits
  resLFC <- lfcShrink(dds, coef=coef, res=res) 
  res$padj_Mine <- p.adjust(res$pvalue, method="BH")
  resLFC$padj_Mine <- p.adjust(resLFC$pvalue, method="BH")
  
  #compute Storey's qvalues
  idx <- which(!is.na(resLFC$pvalue))
  resLFC$qvalue_Mine <- rep(NA,times=length(resLFC$pvalue))
  pvals_temp <- resLFC$pvalue[idx]
  qvals_temp <- qvalue(p = pvals_temp)
  resLFC$qvalue_Mine[idx] <- qvals_temp$qvalues
  
  
  # resLFC <- cbind(Gene, resLFC)
  return(resLFC)
  
}

# ---------- Load config ----------
# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
 stop("Usage: Rscript DESeq2_hybrids.R <config_file.tsv>\n",
      "Example: Rscript DESeq2_hybrids.R DESeq2_hybrids_config.tsv")
}

config_file <- args[1]
cat("Using config file:", config_file, "\n")

config <- read.table(config_file,
                     header = FALSE,
                     sep = "\t",
                     stringsAsFactors = FALSE)

config <- setNames(config$V2, config$V1)

# parameters
input_path     <- config[["input_path"]]
output_path     <- config[["output_path"]]
ref_species    <- config[["ref_species"]]
alt_species    <- config[["alt_species"]]
counts_suffix  <- config[["counts_suffix"]]
metadata_file  <- config[["metadata_file"]]
test_columns   <- strsplit(config[["test_columns"]], ",")[[1]]


# ---------- Load metadata ----------
coldata <- read.table(
  metadata_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = TRUE,
  strip.white = TRUE,
  row.names = "sample"
)

# rownames(coldata) <- coldata$sample

# ---------- Process both species ----------
for (species in c(ref_species, alt_species)) {
  message("Processing species: ", species)
  species_folder <- file.path(input_path, species)
  
  # 1) raw counts for DESeq2
  counts <- load_counts(species_folder, counts_suffix)

  # print head of counts and coldata for debugging
  print("Counts head:")
  print(head(counts))
  print("Coldata head:")
  print(head(coldata))
  
  # 2) TPMs from the same files
  tpms <- load_tpms(species_folder, counts_suffix)
  
  # --- compute mean TPM across samples ---
  ref_tpm_cols   <- grep("_ref_tpm$",   colnames(tpms), value = TRUE)
  alt_tpm_cols   <- grep("_alt_tpm$",   colnames(tpms), value = TRUE)
  total_tpm_cols <- grep("_total_tpm$", colnames(tpms), value = TRUE)
  
  mean_ref_tpm   <- rowMeans(tpms[, ref_tpm_cols,   drop = FALSE])
  mean_alt_tpm   <- rowMeans(tpms[, alt_tpm_cols,   drop = FALSE])
  mean_total_tpm <- rowMeans(tpms[, total_tpm_cols, drop = FALSE])
  # ----------------------------------------
  
  print(colnames(counts))
  print(rownames(coldata))

  # run DESeq2
  dds_res <- run_deseq2(counts, coldata, test_columns)

  output_filename <- paste(output_path, 'DESeq2_aligned2', species, '.txt', sep="")
  output_file <- cbind(dds_res[,"log2FoldChange"],
                       dds_res[,"pvalue"],
                       dds_res[,"padj"],
                       dds_res[,"padj_Mine"],
                       dds_res[,"qvalue_Mine"])
  output_file <- cbind(gene = rownames(output_file), output_file)
  colnames(output_file) <- c("Gene", "log2FoldChange","pvalue","padj","padj_mine","qval_mine")
  
  # ---- add mean TPMs, matched by gene name ----
  # convert matrix → data.frame here
  output_file <- as.data.frame(output_file, stringsAsFactors = FALSE)
  
  output_file$mean_ref_tpm   <- mean_ref_tpm[output_file$Gene]
  output_file$mean_alt_tpm   <- mean_alt_tpm[output_file$Gene]
  output_file$mean_total_tpm <- mean_total_tpm[output_file$Gene]
  
  # (optional) nice column order
  output_file <- output_file[, c("Gene",
                                 "log2FoldChange","pvalue","padj","padj_mine","qval_mine",
                                 "mean_ref_tpm","mean_alt_tpm","mean_total_tpm")]
  # ---------------------------------------------
  
  write.table(output_file,
              file = output_filename,
              sep = "\t",
              quote = FALSE,
              na = "",
              row.names = FALSE,
              col.names = TRUE)
}




