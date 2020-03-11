library(biomaRt)          # useMart(), getBM()
library(dplyr)            # select(), %>%
library(data.table)       # uniqueN()


# #' Search RefSeq mRNA IDs of your interactions data in Ensembl database to obtain matching Ensembl transcript IDs.
# #'
# #' @param interactions Data Frame of miRNA - mRNA interactions. The column name of the mRNAs needs to be "mRNA".
# #' @param ensembl_mart Mart object compatible with the biomaRt query. Gives the possibility to search older versions of the Ensembl database. The default is the current Ensembl database version.
# #' @param id.type Search Ensembl database for "validated", "predicted" or "both" RefSeq mRNA IDs. The default will search for both ID types.
# #' @return interactions Modified interactions data with matching Ensembl transcript IDs
preproc.idconversion.refseqtoensembl <- function(interactions, ensembl_mart = NULL, id.type = 'both'){
  
  
  # Step 1: Prepare reqired data
  refseq_ids <- unique(interactions$mRNA)
  if(is.null(ensembl_mart))
    ensembl_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  
  # Step 2: Retrieve Ensemble IDs via biomaRt
  mappings <- NULL
  additional.attributes <- c('ensembl_transcript_id', 'ensembl_gene_id')
  cols <- c('mRNA', 'ensembl_transcript_id', 'ensembl_gene_id')
  
  if(id.type == 'validated'){
    mappings <- getBM(attributes = c('refseq_mrna', additional.attributes),
                      filters = 'refseq_mrna',
                      values = refseq_ids,
                      mart = ensembl_mart)
    colnames(mappings) <- cols
  } 
  else if(id.type == 'predicted'){
    mappings <- getBM(attributes = c('refseq_mrna_predicted', additional.attributes),
                      filters = 'refseq_mrna_predicted',
                      values = refseq_ids,
                      mart = ensembl_mart)
    colnames(mappings) <- cols
  } 
  else if(id.type == 'both'){
    refseq_ensembl <- getBM(attributes = c('refseq_mrna', additional.attributes),
                            filters = 'refseq_mrna',
                            values = refseq_ids,
                            mart = ensembl_mart)
    predrefseq_ensembl <- getBM(attributes = c('refseq_mrna_predicted', additional.attributes),
                                filters = 'refseq_mrna_predicted',
                                values = refseq_ids,
                                mart = ensembl_mart)
    colnames(refseq_ensembl) <- cols
    colnames(predrefseq_ensembl) <- cols
    mappings <- rbind(refseq_ensembl, predrefseq_ensembl)
  }
  else {
    cat("\nIncorrect value for parameter id_type\nPlease choose between 'validated', 'predicted' and 'both'\nRefSeq to Ensemble IDs conversion stopped...\n")
    return(NULL)
  }
  
  stopifnot(is.null(mappings))
  
  
  # Step 3: Remove entries where no Ensemble IDs were found
  # THIS MIRRETI VERSION USES ONLY DISTINCT MAPPINGS  
  ambiguous_mrnas <- mappings[duplicated(mappings$mRNA), ]
  ambiguous_transcripts <- mappings[duplicated(mappings$ensembl_transcript_id), ]
  distinct_mappings <- mappings[!(mappings$mRNA %in% ambiguous_mrnas$mRNA) &
                                  !(mappings$ensembl_transcript_id %in% ambiguous_transcripts$ensembl_transcript_id), ]

  distinct_interactions <- interactions[interactions$mRNA %in% distinct_mappings$mRNA, ]
  interactions <- merge(distinct_interactions, distinct_mappings, by = c("mRNA"), allow.cartesian = TRUE)
  interactions <- select(interactions, ensembl_gene_id, ensembl_transcript_id, mRNA, miRNA, binding_site, binding_probability, genomic_region)
  
  return(interactions)
    rm(refseq_ids, mappings, refseq_ensembl, predrefseq_ensembl, ambiguous_mrnas, ambiguous_transcripts, distinct_interactions, distinct_mappings, id.type)
}


preproc.filter.samples <- function(tpm_expr, mir_expr, sample_annotation, primary.disease, conditions, log = FALSE){
  start.time <- 0
  end.time <- 0
  
  
  # Step 1: filter samples 
      if(log){
        cat("\n\nMIRRETI preprocess report\n")
        cat("\n##INPUT REPORT##\n")
        cat("\tSamples\n")
        cat(paste("\t\tTotal samples: ", uniqueN(sample_annotation$sample), "\n", sep = ""))
        cat(paste("\t\tNumber of cancer types: ", uniqueN(sample_annotation$X_primary_disease), "\n", sep = ""))
        cat("\tTPM expression data\n")
        cat(paste("\t\tTranscripts: ", uniqueN(rownames(tpm_expr)), "\n", sep = ""))
        cat(paste("\t\tNumbeer of samples: ", uniqueN(colnames(tpm_expr)), "\n", sep = ""))
        cat("\tmiRNA expression data\n")
        cat(paste("\t\tmiRNAs: ", uniqueN(rownames(mir_expr)), "\n", sep = ""))
        cat(paste("\t\tNumber of samples: ", uniqueN(colnames(mir_expr)), "\n", sep = ""))
        cat(paste("\tcancer types: ", paste(primary.disease, collapse = "|"), "\n", sep = ""))
        cat(paste("\tconditions: ", paste(conditions, collapse = "|"), "\n\n", sep = ""))
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart filtering TCGA samples for cancer types and conditions...\n", sep = ""))
      }
  sample_annotation <- sample_annotation %>% filter(X_primary_disease %in% primary.disease &
                                                    sample_type_id %in% conditions)
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\tSUCCESSFUL\n", sep = ""))
        cat(paste("\tTotal samples: ", uniqueN(sample_annotation$sample), "\n", sep = ""))
        cat("!!!!!!!!Hier könnten Ihre Frequenzen stehen!!!!!!!!\n")
        cat(paste("Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
      }
  
  
  # Step 2: subset TPM expression data and preprocess
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart subsetting TPM expression data for samples and filter for NAs and low variance...\n", sep = ""))
      }
  tpm_expr <- tpm_expr[ , colnames(tpm_expr) %in% sample_annotation$sample]
  tpm_expr <- preproc.filter.variance(preproc.filter.na(expr_data = tpm_expr, 
                                                        na.equivalent = min(tpm_expr)))
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\tSUCCESSFUL\n", sep = ""))
        cat(paste("\tTranscripts: ", uniqueN(rownames(tpm_expr)), "\n", sep = ""))
        cat(paste("\tRemaining samples: ", uniqueN(colnames(tpm_expr)), "\n", sep = ""))
        cat(paste("Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
  }
  
  
  # Step 3: subset miRNA expression data and preprocess
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart subsetting miRNA expression data for samples and filter for NAs and low variance...\n", sep = ""))
      }
  mir_expr <- mir_expr[ , colnames(mir_expr) %in% sample_annotation$sample]
  mir_expr <- preproc.filter.variance(preproc.filter.na(expr_data = mir_expr, 
                                                        na.equivalent = min(mir_expr)))
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\tSUCCESSFUL\n", sep = ""))
        cat(paste("\tmiRNAs: ", uniqueN(rownames(mir_expr)), "\n", sep = ""))
        cat(paste("\tRemaining samples: ", uniqueN(colnames(mir_expr)), "\n", sep = ""))
        cat(paste("Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
      }
  
  
  # Step 4: Subsetting sample and expression data for samples included in all datasets
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart subsetting sample and expression data for samples included in all datasets...\n", sep = ""))
      }
  sample_annotation <- sample_annotation %>% filter(sample %in% colnames(tpm_expr) & sample %in% colnames(mir_expr))
  tpm_expr <- tpm_expr[ , colnames(tpm_expr) %in% sample_annotation$sample]
  mir_expr <- mir_expr[ , colnames(mir_expr) %in% sample_annotation$sample]
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\tSUCCESSFUL\n", sep = ""))
        cat(paste("Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
        
        cat("\n##FINAL REPORT##\n")
        cat("\tSamples\n")
        cat(paste("\t\tTotal samples: ", uniqueN(sample_annotation$sample), "\n", sep = ""))
        cat("!!!!!!!!Hier könnten Ihre Frequenzen stehen!!!!!!!!\n")
        cat("TPM expression data\n")
        cat(paste("\t\tTranscripts: ", uniqueN(rownames(tpm_expr)), "\n", sep = ""))
        cat(paste("\t\tRemaining samples: ", uniqueN(colnames(tpm_expr)), "\n", sep = ""))
        cat("miRNA expression data\n")
        cat(paste("\t\tmiRNAs: ", uniqueN(rownames(mir_expr)), "\n", sep = ""))
        cat(paste("\t\tRemaining samples: ", uniqueN(colnames(mir_expr)), "\n\n", sep = ""))
        cat("End of the MIRRETI preprocess report\n\n")
      }
  
  return(list(tpm_brca_01.11, mir_brca_01.11, samples_brca_01.11))
}

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


preproc.filter.actualTranscriptTargets <- function(tpm_expr, mir_expr, interactions){
  mir_expr <- mir_expr[rownames(mir_expr) %in% interactions$miRNA , colnames(mir_expr) %in% colnames(tpm_expr)]
  tpm_expr <- tpm_expr[rownames(tpm_expr) %in% interactions$ensembl_transcript_id, colnames(tpm_expr) %in% colnames(mir_expr)]
  interactions <- interactions[interactions$miRNA %in% rownames(mir_expr) & interactions$ensembl_transcript_id %in% rownames(tpm_expr), ]
  return(list(mir_expr, tpm_expr, interactions))
}