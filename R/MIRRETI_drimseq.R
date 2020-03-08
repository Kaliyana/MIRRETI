library(DRIMSeq)
library(biomaRt)
library(data.table)
library(dplyr)

coredata <- function(){
  tpm.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/tcga_Kallisto_tpm"
  sampleannot.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/TCGA_phenotype_denseDataOnlyDownload.tsv"
  primary.disease <- "breast invasive carcinoma"
  ensembl_mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "http://jan2020.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")
  drimseq <- readRDS("D:/Bioinformatics/Bachelordata/thirdShot/datapacks/drimseq-1.0_breastinvasivecarcinoma_EnsemblGenes99.RDS")
  
  rm(tpm.filepath, sampleannot.filepath)
  rm(primary.disease)
}

mirreti.drimseq <- function(tpm.filepath, sampleannot.filepath, primary.disease, ensembl_mart, log = FALSE){
  start.time <- 0
  end.time <- 0
  
  
  # Step 1: read in data
      if(log){
        start.time <- Sys.time()
        cat("\n\nMIRRETI DRIMSeq report\n")
        cat(paste(start.time, "\tStart reading the data...\n", sep = ""))
      }
  tpm_expr <- read.expression.file(tpm.filepath)
  sample_annotation <- read.sample.annotation.file(sampleannot.filepath)
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\tReading data was successfull\n", sep = ""))
        cat(paste("Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
      }
  
  
  # Step 2: filter samples and expression data
  # samples for cancer type and condition
  # expression data for samples, NAs and variance
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart filtering samples and expression data...\n", sep = ""))
      }
  samples_brca_01.11 <- sample_annotation[sample_annotation$X_primary_disease == primary.disease &
                                            (sample_annotation$sample_type_id == 01 |
                                               sample_annotation$sample_type_id == 11), ]
  samples_brca_01.11 <- select(samples_brca_01.11, sample, sample_type_id)
  colnames(samples_brca_01.11) <- c("sample_id", "condition")
  
  tpm_brca_01.11 <- tpm_expr[ , colnames(tpm_expr) %in% samples_brca_01.11$sample_id]
  tpm_brca_01.11 <- preproc.filter.variance(preproc.filter.na(expr_data = tpm_brca_01.11, 
                                                              na.equivalent = min(tpm_brca_01.11)))
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\tFiltering samples and expression data was successfull\n", sep = ""))
        cat(paste("Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
      }
  
  
  # Step 3: build gene model from Ensembl
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart building gene model on Ensembl transcripts...\n", sep = ""))
      }
  biomart_annot <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"),
                         filters = "ensembl_transcript_id",
                         values = rownames(tpm_brca_01.11),
                         mart = ensembl_mart)
  biomart_annot <- select(biomart_annot, ensembl_gene_id, ensembl_transcript_id)
  colnames(biomart_annot) <- c("gene_id", "feature_id")
  
  biomart_annot_multitranscript <- data.frame(matrix(ncol = ncol(biomart_annot), nrow = 0))
  colnames(biomart_annot_multitranscript) <- colnames(biomart_annot)
  for(g in unique(biomart_annot$gene_id)){
    sub <- biomart_annot[biomart_annot$gene_id == g, ]
    if(uniqueN(sub$feature_id) <= 1)
      next()
    biomart_annot_multitranscript <- rbind(biomart_annot_multitranscript, sub)
  }
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\tBuilding gene model was successfull\n", sep = ""))
        cat(paste("Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
      }
  
  
  # Step 4: prepare samples and counts data for DEA
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart preparing samples and counts data for differential expression...\n", sep = ""))
      }
  tpm_brca_01.11 <- tpm_brca_01.11[rownames(tpm_brca_01.11) %in% biomart_annot_multitranscript$feature_id, ]
  tpm_brca_01.11 <- 2^tpm_brca_01.11
  tpm_brca_01.11$feature_id <- rownames(tpm_brca_01.11)
  counts <- cbind(biomart_annot_multitranscript, tpm_brca_01.11, by = "feature_id")
  rownames(counts) <- NULL
  samples <- samples_brca_01.11[samples_brca_01.11$sample_id %in% colnames(tpm_brca_01.11), ]
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\tPreparing data for DRIMSeq was successfull\n", sep = ""))
        cat(paste("Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
      }
  
  # Step 5: DRIMSeq
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart running DRIMSeq...\n", sep = ""))
      }
  drimseq <- dmDSdata(counts=counts, samples=samples)
  #table(table(counts(d)$gene_id))
  design_full <- model.matrix(~condition, data=DRIMSeq::samples(drimseq))
  #colnames(design_full)
  set.seed(01101990)
  system.time({
    drimseq <- dmPrecision(drimseq, design=design_full)
    drimseq <- dmFit(drimseq, design=design_full)
    drimseq <- dmTest(drimseq, coef="condition")
  })
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\tRunning DRIMSeq was successfull\n", sep = ""))
        cat(paste("Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
        cat("End of the MIRRETI DRIMSeq report")
      }
  
  
  return(drimseq)
}

mirreti.drimseq.results <- function(drimseq){
  plotPValues(drimseq)
  drimseq@samples[, "condition"] <- factor(drimseq@samples[, "condition"])
  
  results <- DRIMSeq::results(drimseq)
  results <- results[order(results$pvalue, decreasing = FALSE), ]
  results_txp <- DRIMSeq::results(drimseq, level = "feature")
  no.na <- function(x) ifelse(is.na(x), 1, x)
  results$pvalue <- no.na(results$pvalue)
  results_txp$pvalue <- no.na(results_txp$pvalue)
  
  idx <- which(results$adj_pvalue < 0.05)[1]
  plotProportions(drimseq, 
                  gene_id = "ENSG00000185624", 
                  group_variable = "condition")
}


read.expression.file <- function(expression.filepath){
  expr_data <- data.frame(fread(expression.filepath, header = T, sep = "\t"), row.names = 1, check.names = F)
  rownames(expr_data) <- lapply(strsplit(rownames(expr_data), "\\."), '[', 1)
  return(expr_data)
}

read.sample.annotation.file <- function(sampleannotation.filepath){
  sample_annotation <- read.csv(sampleannotation.filepath, header = T, sep = "\t")
  return(sample_annotation)
}

# na.equivalent transcripts:  -9.9658
preproc.filter.na <- function(expr_data, na.equivalent, max.na.percent = 0.2){
  na_content <- rowSums(expr_data == na.equivalent)
  expr_data <- expr_data[na_content/ncol(expr_data) <= max.na.percent, ]
  na_content <- colSums(expr_data == na.equivalent)
  expr_data <- expr_data[ , na_content/nrow(expr_data) <= max.na.percent]
  return(expr_data)
}

preproc.filter.variance <- function(expr_data, var.threshold = 0.2){
  expr_data$variance <- apply(data.frame(expr_data), 1, var)
  expr_data <- expr_data[order(expr_data$variance, decreasing = T), ]
  threshold <- floor((1-var.threshold) * nrow(expr_data))
  expr_data <- expr_data[1:threshold, colnames(expr_data) != 'variance']
  return(expr_data)
}

mirreti.drimseq.results(drimseq)

