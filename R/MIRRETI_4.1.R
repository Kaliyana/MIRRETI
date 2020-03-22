#load("~/bachelor/data/core_data/mirreti_coredata.RDS")
set.seed(01101990)
library(data.table)                             # uniqueN(), unique()
library(dplyr)                                  # %>%, select()
library(biomaRt)                                # useMart(), getBM(), listMarts()
library(doParallel)                             # makePSOCKcluster(), registerDoParallel(), stopCluster()
library(SPONGE)                                 # sponge_gene_miRNA_interaction_filter()
library(DRIMSeq)                                # dmDSdata(), model.matrix(), DRIMSeq::samples(), dmPrecision(), dmFit(), dmTest()



#------------------------------------------------------------------------------------------------
#                           Flexible subsetting of TCGA expression data
#------------------------------------------------------------------------------------------------

# #' DESCRIBE METHODE
# #'
# #' @param tpm_expr
# #' @param mir_expr
# #' @param sample_annotation
# #' @param primary.diseases
# #' @param sample.type.ids
# #' @param max.na.percent
# #' @param log
# #' @return tpm_expr, mir_expr, sample_annotation
tcgaSamples.subset <- function(tpm_expr, mir_expr, sample_annotation, primary.diseases, sample.type.ids = "1", max.na.percent = 0.2, log = FALSE, log.file = NULL){
      start.time <- 0
      end.time <- 0
      if(!is.null(log.file)){
        sink(log.file, append = TRUE)
      }
      
  # submethods
  exprdata.filter.na <- function(expr_data, max.na.percent = 0.2){
        na_content <- rowSums(expr_data == min(expr_data))
        expr_data <- expr_data[na_content/ncol(expr_data) <= max.na.percent, ]
        na_content <- colSums(expr_data == min(expr_data))
        expr_data <- expr_data[ , na_content/nrow(expr_data) <= max.na.percent]
        return(expr_data)
      }
  
  # Step 1: filter samples 
      if(log){
        cat("______________________________________________________________________\n\n")
        cat(paste("\U1F183\U1F172\U1F176\U1F170", "EXPRESSION DATA SUBSET REPORT\n"))
        cat("\n  \U1F664 \U1F158\U1F15D\U1F15F\U1F164\U1F163 \U1F161\U1F154\U1F15F\U1F15E\U1F161\U1F163 \U1F666\n")
        cat("  SAMPLES\n")
        cat(paste("    Total samples: ", uniqueN(sample_annotation$sample), "\n", sep = ""))
        cat(paste("    number of cancer types: ", uniqueN(sample_annotation$X_primary_disease), "\n", sep = ""))
        cat("  TRANSCRIPT EXPRESSION DATA\n")
        cat(paste("    transcripts: ", uniqueN(rownames(tpm_expr)), "\n", sep = ""))
        cat(paste("    number of samples: ", uniqueN(colnames(tpm_expr)), "\n", sep = ""))
        cat("  MICRORNA EXPRESSION DATA\n")
        cat(paste("    microRNAs: ", uniqueN(rownames(mir_expr)), "\n", sep = ""))
        cat(paste("    number of samples: ", uniqueN(colnames(mir_expr)), "\n", sep = ""))
        cat(paste("  CANCER TYPES: ", paste(primary.diseases, collapse = "|"), "\n", sep = ""))
        cat(paste("  SAMPLE TYPES: ", paste(sample.type.ids, collapse = "|"), "\n\n\n", sep = ""))
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart filtering TCGA samples for cancer and sample types...\n", sep = ""))
      }
  sample_annotation <- sample_annotation %>% filter(X_primary_disease %in% primary.diseases &
                                                      sample_type_id %in% sample.type.ids)
  sample_annotation$condition <- paste(gsub(" ", "", sample_annotation$X_primary_disease), sample_annotation$sample_type_id, sep = "_")
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\t\U2705 \U1F1F8\U1F1FA\U1F1E8\U1F1E8\U1F1EA\U1F1F8\U1F1F8\U1F1EB\U1F1FA\U1F1F1\n", sep = ""))
        cat(paste("\tTotal samples: ", uniqueN(sample_annotation$sample), "\n", sep = ""))
        for(d in unique(sample_annotation$X_primary_disease)){
          d_sub <- sample_annotation[sample_annotation$X_primary_disease == d, ]
          cat(paste("\t", d, ": ", nrow(d_sub), "\n", sep = ""))
          for(st in unique(sample_annotation$sample_type_id)){
            cat(paste("\t - sample type ", st, ": ", nrow(d_sub[d_sub$sample_type_id == st, ]), "\n", sep = ""))
          }
        }
        cat(paste("\U231B Time required:\t", round(difftime(end.time, start.time, units = "secs"), digits = 3), " secs\n\n", sep = ""))
      }
  
  
  # Step 2: subset TPM expression data and filter transcripts and samples for NAs
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart subsetting transcript expression data for samples and filter for NAs...\n", sep = ""))
      }
  tpm_expr <- tpm_expr[ , colnames(tpm_expr) %in% sample_annotation$sample]
  tpm_expr <- exprdata.filter.na(tpm_expr, max.na.percent)
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\t\U2705 \U1F1F8\U1F1FA\U1F1E8\U1F1E8\U1F1EA\U1F1F8\U1F1F8\U1F1EB\U1F1FA\U1F1F1\n", sep = ""))
        cat(paste("\ttranscripts: ", uniqueN(rownames(tpm_expr)), "\n", sep = ""))
        cat(paste("\tremaining samples: ", uniqueN(colnames(tpm_expr)), "\n", sep = ""))
        cat(paste("\U231B Time required:\t", round(difftime(end.time, start.time, units = "secs"), digits = 3), " secs\n\n", sep = ""))
      }
  
  
  # Step 3: subset miRNA expression data and filter microRNAs and samples for NAs
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart subsetting microRNA expression data for samples and filter for NAs...\n", sep = ""))
      }
  mir_expr <- mir_expr[ , colnames(mir_expr) %in% sample_annotation$sample]
  mir_expr <- exprdata.filter.na(mir_expr, max.na.percent)
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\t\U2705 \U1F1F8\U1F1FA\U1F1E8\U1F1E8\U1F1EA\U1F1F8\U1F1F8\U1F1EB\U1F1FA\U1F1F1\n", sep = ""))
        cat(paste("\tmicroRNAs: ", uniqueN(rownames(mir_expr)), "\n", sep = ""))
        cat(paste("\tremaining samples: ", uniqueN(colnames(mir_expr)), "\n", sep = ""))
        cat(paste("\U231B Time required:\t", round(difftime(end.time, start.time, units = "secs"), digits = 3), " secs\n\n", sep = ""))
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
        cat(paste(end.time, "\t\U2705 \U1F1F8\U1F1FA\U1F1E8\U1F1E8\U1F1EA\U1F1F8\U1F1F8\U1F1EB\U1F1FA\U1F1F1\n", sep = ""))
        cat(paste("\U231B Time required:\t", round(difftime(end.time, start.time, units = "secs"), digits = 3), " secs\n\n", sep = ""))
        cat("\n  \U1F664 \U1F15E\U1F164\U1F163\U1F15F\U1F164\U1F163 \U1F161\U1F154\U1F15F\U1F15E\U1F161\U1F163 \U1F666\n")
        cat("  SAMPLES\n")
        cat(paste("    Total samples: ", uniqueN(sample_annotation$sample), "\n", sep = ""))
        for(d in unique(sample_annotation$X_primary_disease)){
          d_sub <- sample_annotation[sample_annotation$X_primary_disease == d, ]
          cat(paste("    ", d, ": ", nrow(d_sub), "\n", sep = ""))
          for(st in unique(sample_annotation$sample_type_id)){
            cat(paste("     - sample type ", st, ": ", nrow(d_sub[d_sub$sample_type_id == st, ]), "\n", sep = ""))
          }
        }
        cat("  TRANSCRIPT EXPRESSION DATA\n")
        cat(paste("    transcripts: ", uniqueN(rownames(tpm_expr)), "\n", sep = ""))
        cat(paste("    remaining samples: ", uniqueN(colnames(tpm_expr)), "\n", sep = ""))
        cat("  MICRORNA EXPRESSION DATA\n")
        cat(paste("    microRNAs: ", uniqueN(rownames(mir_expr)), "\n", sep = ""))
        cat(paste("    remaining samples: ", uniqueN(colnames(mir_expr)), "\n\n", sep = ""))
        cat("END OF THE TCGA EXPRESSION DATA SUBSET REPORT\n")
        cat("______________________________________________________________________\n")
      }
  
  sample_annotation <- dplyr::select(sample_annotation, sample, condition)
  closeAllConnections()
  return(list(tpm_expr, mir_expr, sample_annotation))
}

# no submethod
tcgaSamples.analyse.sampleAnnotation <- function(sample_annotation){
  table(sample_annotation$X_primary_disease)[order(table(sample_annotation$X_primary_disease), decreasing = TRUE)]
}



#------------------------------------------------------------------------------------------------
#             miRWalk specific ID conversion from RefSeq mRNA to Ensembl trancripts
#------------------------------------------------------------------------------------------------

# #' Search RefSeq mRNA IDs of your miRWalk interactions data in Ensembl database to obtain matching Ensembl transcript IDs.
# #'
# #' @param interactions: Data Frame of miRNA - mRNA interactions. The column name of the mRNAs needs to be "mRNA".
# #' @param ensembl_mart: Mart object compatible with the biomaRt query. Gives the possibility to search older versions of the Ensembl database. The default is the current Ensembl database version.
# #' @param id.type: Search Ensembl database for "validated", "predicted" or "both" RefSeq mRNA IDs. The default will search for both ID types.
# #' @return interactions: Modified interactions data with matching Ensembl transcript IDs
mirwalkInteractions.idconversion <- function(interactions, ensembl_mart = NULL, id.type = 'both', log = FALSE){
      start.time <- 0
      end.time <- 0
  if(is.null(ensembl_mart)){
    ensembl_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  }
      
  # Step 1: Retrieve Ensemble IDs via biomaRt
      if(log){
        cat("______________________________________________________________________\n\n")
        cat(paste("\U1F17C\U1F178\U1F181\U1F186\U1F170\U1F17B\U1F17A", "INTERACTIONS ID CONVERSION REPORT\n"))
        cat("\n  \U1F664 \U1F158\U1F15D\U1F15F\U1F164\U1F163 \U1F161\U1F154\U1F15F\U1F15E\U1F161\U1F163 \U1F666\n")
        cat("  INTERACTIONS\n")
        cat(paste("    Total interactions: ", nrow(interactions), "\n", sep = ""))
        cat(paste("    microRNAs: ", uniqueN(interactions$miRNA), "\n", sep = ""))
        cat(paste("    RefSeq mRNAs: ", uniqueN(interactions$mRNA), "\n", sep = ""))
        cat(paste("    genes: ", uniqueN(interactions$Genesymbol), "\n", sep = ""))
        cat(paste("  ENSEMBL MART: ", listMarts(ensembl_mart)[1,2], "\n", sep = ""))
        cat(paste("  ID TYPE: ", id.type, "\n\n\n", sep = ""))
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart retreiving Ensembl annotation for RefSeq mRNA IDs...\n", sep = ""))
      }
  refseq_ids <- unique(interactions$mRNA)
  mappings <- NULL
  additional.attributes <- c('ensembl_transcript_id', 'ensembl_gene_id')
  cols <- c('mRNA', additional.attributes)
  
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
      rm(refseq_ensembl, predrefseq_ensembl)
      gc()
  }
  else {
    cat("\nIncorrect value for parameter id_type\nPlease choose between 'validated', 'predicted' and 'both'\nRefSeq to Ensemble IDs conversion stopped...\n")
    stopifnot(FALSE)
  }
  stopifnot(!is.null(mappings))
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\t\U2705 \U1F1F8\U1F1FA\U1F1E8\U1F1E8\U1F1EA\U1F1F8\U1F1F8\U1F1EB\U1F1FA\U1F1F1\n", sep = ""))
        cat(paste("\U231B Time required:\t", round(difftime(end.time, start.time, units = "secs"), digits = 3), " secs\n\n", sep = ""))
      }
  
  
  # Step 2: Analysis of the found ID mappings
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart analysing found RefSeq to Ensembl ID mappings...\n", sep = ""))
      }
  not_mapped_interactions <- interactions[!interactions$mRNA %in% mappings$mRNA, ]
  mapped_interactions <- interactions[interactions$mRNA %in% mappings$mRNA, ]
  
  ambiguous_mrnas <- mappings[duplicated(mappings$mRNA), 'mRNA']
  ambiguous_transcripts <- mappings[duplicated(mappings$ensembl_transcript_id), 'ensembl_transcript_id']
  ambiguous_mappings <- mappings[mappings$mRNA %in% ambiguous_mrnas | mappings$ensembl_transcript_id %in% ambiguous_transcripts, ]
  ambiguous_interactions <- interactions[interactions$mRNA %in% ambiguous_mappings$mRNA, ]
  
  distinct_mappings <- mappings[!mappings$mRNA %in% ambiguous_mrnas & !mappings$ensembl_transcript_id %in% ambiguous_transcripts, ]
  distinct_interactions <- interactions[interactions$mRNA %in% distinct_mappings$mRNA, ]
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\t\U2705 \U1F1F8\U1F1FA\U1F1E8\U1F1E8\U1F1EA\U1F1F8\U1F1F8\U1F1EB\U1F1FA\U1F1F1\n", sep = ""))
        cat(paste("\U231B Time required:\t", round(difftime(end.time, start.time, units = "secs"), digits = 3), " secs\n\n", sep = ""))
        cat("  \U1F664 \U1F15C\U1F150\U1F15F\U1F15F\U1F158\U1F15D\U1F156 \U1F161\U1F154\U1F15F\U1F15E\U1F161\U1F163 \U1F666\n")
        cat("  NOT MAPPED INTERACTIONS\n")
        cat(paste("    Total not mapped interactions: ", nrow(not_mapped_interactions), "\n", sep = ""))
        cat(paste("    mRNAs: ", uniqueN(not_mapped_interactions$mRNA)), "\n", sep = "")
        cat("  AMBIGUOUSLY MAPPED INTERACTIONS\n")
        cat(paste("    Total ambiguously mapped interactions: ", nrow(ambiguous_interactions), "\n", sep = ""))
        cat(paste("    mRNAs: ", uniqueN(ambiguous_interactions$mRNA)), "\n", sep = "")
        cat("  DISTINCTLY MAPPED INTERACTIONS\n")
        cat(paste("    Total distinctly mapped interactions: ", nrow(distinct_interactions), "\n", sep = ""))
        cat(paste("    mRNAs: ", uniqueN(distinct_interactions$mRNA), "\n", sep = ""))
        cat(paste("    genes: ", uniqueN(distinct_interactions$Genesymbol), "\n", sep = ""))
        cat(paste("     -> genes with removed transcripts: ", uniqueN(distinct_interactions[distinct_interactions$Genesymbol %in% not_mapped_interactions$Genesymbol | distinct_interactions$Genesymbol %in% ambiguous_interactions$Genesymbol, 'Genesymbol']), "\n", sep = ""))
        cat(paste("     -> totally removed genes: ", uniqueN(interactions[!interactions$Genesymbol %in% distinct_interactions$Genesymbol, 'Genesymbol']), "\n\n", sep = ""))
      }
    rm(mappings, not_mapped_interactions, mapped_interactions, ambiguous_mrnas, ambiguous_transcripts, ambiguous_mappings, ambiguous_interactions)
    gc()
  
    
  # Step 3: Combine annotations of distinctly mapped IDs
  # THIS MIRRETI VERSION USES ONLY DISTINCT MAPPINGS
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart combining annotations of distinctly mapped IDs...\n", sep = ""))
      }
  interactions <- merge(distinct_interactions, distinct_mappings, by = c("mRNA"), allow.cartesian = TRUE)
  interactions <- dplyr::select(interactions, all_of(additional.attributes), mRNA, miRNA, binding_site, binding_probability, genomic_region)
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\t\U2705 \U1F1F8\U1F1FA\U1F1E8\U1F1E8\U1F1EA\U1F1F8\U1F1F8\U1F1EB\U1F1FA\U1F1F1\n", sep = ""))
        cat(paste("\U231B Time required:\t", round(difftime(end.time, start.time, units = "secs"), digits = 3), " secs\n\n", sep = ""))
        cat("\n  \U1F664 \U1F15E\U1F164\U1F163\U1F15F\U1F164\U1F163 \U1F161\U1F154\U1F15F\U1F15E\U1F161\U1F163 \U1F666\n")
        cat("  INTERACTIONS\n")
        cat(paste("    Total interactions: ", nrow(interactions), "\n", sep = ""))
        cat(paste("    microRNAs: ", uniqueN(interactions$miRNA), "\n", sep = ""))
        cat(paste("    RefSeq mRNAs: ", uniqueN(interactions$mRNA), "\n", sep = ""))
        cat(paste("    Ensembl transcripts: ", uniqueN(interactions$ensembl_transcript_id), "\n", sep = ""))
        cat(paste("    Ensembl genes: ", uniqueN(interactions$ensembl_gene_id), "\n\n", sep = ""))
        cat("END OF THE MIRWALK INTERACTIONS ID CONVERSION REPORT\n")
        cat("______________________________________________________________________\n")
      }
  
  return(interactions)
}



#------------------------------------------------------------------------------------------------
#                               Crossfilter to align data dimensions
#------------------------------------------------------------------------------------------------

# #' DESCRIBE METHODE
# #'
# #' @param tpm_expr
# #' @param mir_expr
# #' @param interactions
# #' @param var.threshold
# #' @return tpm_expr, mir_expr, interactions
mirretiData.filter.actualTargets <- function(tpm_expr, mir_expr, interactions, var.threshold = NULL, log = FALSE){
      start.time <- 0
      end.time <- 0
      
  # submethods
  exprdata.filter.variance <- function(expr_data, var.threshold = 0.2){
        expr_data$variance <- apply(data.frame(expr_data), 1, var)
        expr_data <- expr_data[order(expr_data$variance, decreasing = T), ]
        threshold <- floor((1-var.threshold) * nrow(expr_data))
        expr_data <- expr_data[1:threshold, colnames(expr_data) != 'variance']
        return(expr_data)
      }
  
  # Step 1: filter out var.threashold least variant columns (molecules) in expression data
      if(log){
        cat("______________________________________________________________________\n\n")
        cat(paste("\U1F17C\U1F178\U1F181\U1F181\U1F174\U1F183\U1F178", "ACTUAL DATA TARGETS REPORT\n"))
        cat("\n  \U1F664 \U1F158\U1F15D\U1F15F\U1F164\U1F163 \U1F161\U1F154\U1F15F\U1F15E\U1F161\U1F163 \U1F666\n")
        cat("  INTERACTIONS\n")
        cat(paste("    Total interactions: ", nrow(interactions), "\n", sep = ""))
        cat(paste("    transcripts: ", uniqueN(interactions$ensembl_transcript_id), "\n", sep = ""))
        cat(paste("    genes: ", uniqueN(interactions$ensembl_gene_id), "\n", sep = ""))
        cat(paste("    microRNAs: ", uniqueN(interactions$miRNA), "\n", sep = ""))
        cat("  TRANSCRIPT EXPRESSION DATA\n")
        cat(paste("    transcripts: ", uniqueN(rownames(tpm_expr)), "\n", sep = ""))
        cat(paste("    number of samples: ", uniqueN(colnames(tpm_expr)), "\n", sep = ""))
        cat("  MICRORNA EXPRESSION DATA\n")
        cat(paste("    microRNAs: ", uniqueN(rownames(mir_expr)), "\n", sep = ""))
        cat(paste("    number of samples: ", uniqueN(colnames(mir_expr)), "\n", sep = ""))
        cat(paste("  VARIANCE THRESHOLD: ", var.threshold, "\n\n\n", sep = ""))
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart filtering expression data for least variant transcript and microRNAs...\n", sep = ""))
      }
  rownames(tpm_expr) <- lapply(strsplit(rownames(tpm_expr), "\\."), '[', 1)
  if(!is.null(var.threshold)){
    tpm_expr <- exprdata.filter.variance(expr_data = tpm_expr,
                                         var.threshold = var.threshold)
    mir_expr <- exprdata.filter.variance(expr_data = mir_expr,
                                         var.threshold = var.threshold)
  }
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\t\U2705 \U1F1F8\U1F1FA\U1F1E8\U1F1E8\U1F1EA\U1F1F8\U1F1F8\U1F1EB\U1F1FA\U1F1F1\n", sep = ""))
        cat(paste("\ttranscripts: ", uniqueN(rownames(tpm_expr)), "\n", sep = ""))
        cat(paste("\tmicroRNAs: ", uniqueN(rownames(mir_expr)), "\n", sep = ""))
        cat(paste("\U231B Time required:\t", round(difftime(end.time, start.time, units = "secs"), digits = 3), " secs\n\n", sep = ""))
      }
  
  # Step 2: match actual overlapping microRNAs and transcript in expression and interactions data
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart subsetting interactions and expression data for commonly included transcripts and microRNAs...\n", sep = ""))
      }
  tpm_expr <- tpm_expr[rownames(tpm_expr) %in% interactions$ensembl_transcript_id, colnames(tpm_expr) %in% colnames(mir_expr)]
  mir_expr <- mir_expr[rownames(mir_expr) %in% interactions$miRNA , colnames(mir_expr) %in% colnames(tpm_expr)]
  interactions <- interactions[interactions$miRNA %in% rownames(mir_expr) & interactions$ensembl_transcript_id %in% rownames(tpm_expr), ]
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\t\U2705 \U1F1F8\U1F1FA\U1F1E8\U1F1E8\U1F1EA\U1F1F8\U1F1F8\U1F1EB\U1F1FA\U1F1F1\n", sep = ""))
        cat(paste("\U231B Time required:\t", round(difftime(end.time, start.time, units = "secs"), digits = 3), " secs\n\n", sep = ""))
        cat("\n  \U1F664 \U1F15E\U1F164\U1F163\U1F15F\U1F164\U1F163 \U1F161\U1F154\U1F15F\U1F15E\U1F161\U1F163 \U1F666\n")
        cat("  INTERACTIONS\n")
        cat(paste("    Total interactions: ", nrow(interactions), "\n", sep = ""))
        cat(paste("    transcripts: ", uniqueN(interactions$ensembl_transcript_id), "\n", sep = ""))
        cat(paste("    genes: ", uniqueN(interactions$ensembl_gene_id), "\n", sep = ""))
        cat(paste("    microRNAs: ", uniqueN(interactions$miRNA), "\n", sep = ""))
        cat("  TRANSCRIPT EXPRESSION DATA\n")
        cat(paste("    transcripts: ", uniqueN(rownames(tpm_expr)), "\n", sep = ""))
        cat(paste("    number of samples: ", uniqueN(colnames(tpm_expr)), "\n", sep = ""))
        cat("  MICRORNA EXPRESSION DATA\n")
        cat(paste("    microRNAs: ", uniqueN(rownames(mir_expr)), "\n", sep = ""))
        cat(paste("    number of samples: ", uniqueN(colnames(mir_expr)), "\n\n", sep = ""))
        cat("END OF THE ACTUAL DATA TARGETS REPORT\n")
        cat("______________________________________________________________________\n")
      }
  
  return(list(tpm_expr, mir_expr, interactions))
}



# #' DESCRIBE METHODE
# #'
# #' @param tpm_expr
# #' @param mir_expr
# #' @param interactions
# #' @param f.test
# #' @param cluster.size
# #' @param log
# #' @return interactions
mirreti.oneCondition <- function(tpm_expr, mir_expr, interactions, f.test = FALSE, cluster.size, log = FALSE){
      cat("\n\n\t\t\t\t*** Start MIRRETI 4.1 one condition analysis ***\n")

  # Step 1: run SPONGE interaction filter and determin BCTs (Binding site Controlled Transcripts)
  interactions <- sponge.correlateData(tpm_expr = tpm_expr,
                                       mir_expr = mir_expr,
                                       interactions = interactions,
                                       cluster.size = cluster.size,
                                       log = log)
  
  # Step 2: determin BCMTPs (Binding site Controlled and Missing Transcript Pairs)
  bcmtp <- mirreti.determin.bcmtp(interactions = interactions,
                                  log = log)
  data <- mirretiData.filter.actualTargets(tpm_expr = tpm_expr, 
                                           mir_expr = mir_expr, 
                                           interactions = bcmtp,
                                           log = log)
  
  # Step 3: correlate expression data with BMTs
  if(f.test %in% c(TRUE, FALSE)){
    bcmtp <- sponge.correlateData(tpm_expr = data[[1]],
                                  mir_expr = data[[2]],
                                  interactions = data[[3]],
                                  f.test = f.test,
                                  cluster.size = ceiling(cluster.size/10),
                                  log = log)
    mirreti.plot.bcmtpStatistics(bcmtp = bcmtp,
                                 log = log)
  }
  else if(f.test == "both"){
    bcmtp_noftest <- sponge.correlateData(tpm_expr = data[[1]],
                                          mir_expr = data[[2]],
                                          interactions = data[[3]],
                                          f.test = FALSE,
                                          cluster.size = ceiling(cluster.size/10),
                                          log = log)
    mirreti.plot.bcmtpStatistics(bcmtp = bcmtp_noftest,
                                 log = log)
    
    bcmtp_ftest <- sponge.correlateData(tpm_expr = data[[1]],
                                        mir_expr = data[[2]],
                                        interactions = data[[3]],
                                        f.test = TRUE,
                                        cluster.size = ceiling(cluster.size/10),
                                        log = log)
    mirreti.plot.bcmtpStatistics(bcmtp = bcmtp_ftest,
                                 log = log)
  }

  return(interactions)
}



# submethod
# anker for variation
mirreti.determin.bcmtp <- function(interactions, log = FALSE){
      start.time <- 0
      end.time <- 0
  
      if(log){
        cat("______________________________________________________________________\n\n")
        cat(paste("\U1F17C\U1F178\U1F181\U1F181\U1F174\U1F183\U1F178", "BINDING SITE CONTROLLED OR MISSING TRANSCRIPT REPORT\n"))
        cat("\n  \U1F664 \U1F158\U1F15D\U1F15F\U1F164\U1F163 \U1F161\U1F154\U1F15F\U1F15E\U1F161\U1F163 \U1F666\n")
        cat("  INTERACTIONS\n")
        cat(paste("    Total interactions: ", nrow(interactions), "\n", sep = ""))
        cat(paste("     - correlated: ", nrow(interactions[interactions$correlation != '-', ]), "\n", sep = ""))
        cat(paste("    transcripts: ", uniqueN(interactions$ensembl_transcript_id), "\n", sep = ""))
        cat(paste("     - correlated: ", uniqueN(interactions[interactions$correlation != '-', 'ensembl_transcript_id']), "\n", sep = ""))
        cat(paste("    genes: ", uniqueN(interactions$ensembl_gene_id), "\n", sep = ""))
        cat(paste("     - correlated: ", uniqueN(interactions[interactions$correlation != '-', 'ensembl_gene_id']), "\n", sep = ""))
        cat(paste("    microRNAs: ", uniqueN(interactions$miRNA), "\n", sep = ""))
        cat(paste("     - correlated: ", uniqueN(interactions[interactions$correlation != '-', 'miRNA']), "\n\n\n", sep = ""))
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart filtering expression data for least variant transcript and microRNAs...\n", sep = ""))
      }
  bcmtp <- data.table(matrix(ncol = 4, nrow = 0))
  colnames(bcmtp) <- c("ensembl_gene_id", "ensembl_transcript_id", "miRNA", "binding_site_event")
    count <- 0
    pb <- txtProgressBar(min = 0, max = uniqueN(interactions$ensembl_gene_id), style = 3)
  for(g in unique(interactions$ensembl_gene_id)){
    g_sub <- interactions[interactions$ensembl_gene_id == g, ]
    t.numb <- uniqueN(g_sub$ensembl_transcript_id)
    if(t.numb < 2){
      next
    }
    
    match <- g_sub[g_sub$correlation != "-", ]
    match <- unique(match$miRNA)
    if(length(match) == 0){
      next
    }
    for(mi in match){
      mi_sub <- g_sub[g_sub$miRNA == mi, ]
      if(uniqueN(mi_sub$ensembl_transcript_id) == t.numb){
        next
      }
      
      bs_losing_transcripts <- unique(g_sub$ensembl_transcript_id[!g_sub$ensembl_transcript_id %in% unique(mi_sub$ensembl_transcript_id)])
      df <- data.table(ensembl_gene_id = rep(g, length(bs_losing_transcripts)),
                       ensembl_transcript_id = bs_losing_transcripts,
                       miRNA = rep(mi, length(bs_losing_transcripts)),
                       binding_site_event = rep("BMT", length(bs_losing_transcripts)))
      bcmtp <- rbind(bcmtp, df)
      
      mi_sub <- mi_sub[mi_sub$correlation != "-", ]
      df <- data.table(ensembl_gene_id = rep(g, uniqueN(mi_sub$ensembl_transcript_id)),
                       ensembl_transcript_id = unique(mi_sub$ensembl_transcript_id),
                       miRNA = rep(mi, uniqueN(mi_sub$ensembl_transcript_id)),
                       binding_site_event = rep("BCT", uniqueN(mi_sub$ensembl_transcript_id)))
      bcmtp <- rbind(bcmtp, df)
    }
    
    count <- count + 1
    setTxtProgressBar(pb, count)
  }
    close(pb)
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\t\U2705 \U1F1F8\U1F1FA\U1F1E8\U1F1E8\U1F1EA\U1F1F8\U1F1F8\U1F1EB\U1F1FA\U1F1F1\n", sep = ""))
        cat(paste("\U231B Time required:\t", round(difftime(end.time, start.time, units = "secs"), digits = 3), " secs\n\n", sep = ""))
        cat("\n  \U1F664 \U1F15E\U1F164\U1F163\U1F15F\U1F164\U1F163 \U1F161\U1F154\U1F15F\U1F15E\U1F161\U1F163 \U1F666\n")
        cat("  BINDING SITE CONTROLLED AND MISSING TRANSCRIPT PAIRS\n")
        cat(paste("    Total pairs: ", nrow(bcmtp), "\n", sep = ""))
        cat(paste("     - bct: ", nrow(bcmtp[bcmtp$binding_site_event == "BCT", ]), "\n", sep = ""))
        cat(paste("     - bmt: ", nrow(bcmtp[bcmtp$binding_site_event == "BMT", ]), "\n", sep = ""))
        cat(paste("    transcripts: ", uniqueN(bcmtp$ensembl_transcript_id), "\n", sep = ""))
        cat(paste("     - bct: ", uniqueN(bcmtp[bcmtp$binding_site_event == "BCT", 'ensembl_transcript_id']), "\n", sep = ""))
        cat(paste("     - bmt: ", uniqueN(bcmtp[bcmtp$binding_site_event == "BMT", 'ensembl_transcript_id']), "\n", sep = ""))
        cat("     \U270E transcripts may be BCT for one microRNA but BMT for another microRNA\n")
        cat(paste("    microRNAs: ", uniqueN(bcmtp$miRNA), "\n", sep = ""))
        cat(paste("    genes: ", uniqueN(bcmtp$ensembl_gene_id), "\n\n", sep = ""))
        cat("END OF THE BINDING SITE CONTROLLED OR MISSING TRANSCRIPT REPORT\n")
        cat("______________________________________________________________________\n")
      }
  
  return(bcmtp)
    rm(df, g_sub, mi_sub, pb, bs_losing_transcripts, count, g, match, mi, t.numb)
}

# submethod
mirreti.plot.bcmtpStatistics <- function(bcmtp, log = FALSE){
  
  # submethods
  determin.geneMeans <- function(bcmtp){
    gene_correlation <- data.table(matrix(nrow = 0, ncol = 1))
    if('coefficient' %in% colnames(bcmtp)){
      gene_correlation <- data.table(matrix(nrow = 0, ncol = 4))
      colnames(gene_correlation) <- c("ensembl_gene_id", "miRNA", "bct_mean_coefficient", "bmt_mean_coefficient")
    }
    else if('fstats' %in% colnames(bcmtp) & 'pval' %in% colnames(bcmtp) & 'p.adj' %in% colnames(bcmtp)){
      gene_correlation <- data.table(matrix(nrow = 0, ncol = 8))
      colnames(gene_correlation) <- c("ensembl_gene_id", "miRNA", "bct_mean_fstats", "bmt_mean_fstats", "bct_mean_pval", "bmt_mean_pval", "bct_mean_p.adj", "bmt_mean_p.adj")
    }
    
    for(g in unique(bcmtp$ensembl_gene_id)){
      g_sub <- bcmtp[bcmtp$ensembl_gene_id == g, ]
      for(mi in unique(g_sub$miRNA)){
        mi_sub <- g_sub[g_sub$miRNA == mi, ]
        if('coefficient' %in% colnames(mi_sub) & all(mi_sub$coefficient == 0)){
          next
        }
        else if('fstats' %in% colnames(bcmtp) & 'pval' %in% colnames(bcmtp) & 'p.adj' %in% colnames(bcmtp) &
                all(mi_sub$fstats == 0) & all(mi_sub$pval == 0) & all(mi_sub$p.adj == 0)){
          next
        }
        
        bct <- mi_sub[mi_sub$binding_site_event == "BCT", ]
        bmt <- mi_sub[mi_sub$binding_site_event == "BMT", ]
        if('coefficient' %in% colnames(mi_sub)){
          dt <- data.table(ensembl_gene_id = g,
                           miRNA = mi,
                           bct_mean_coefficient = mean(bct$coefficient),
                           bmt_mean_coefficient = mean(bmt$coefficient))
          gene_correlation <- rbind(gene_correlation, dt)
        }
        else if('fstats' %in% colnames(bcmtp) & 'pval' %in% colnames(bcmtp) & 'p.adj' %in% colnames(bcmtp)){
          dt <- data.table(ensembl_gene_id = g,
                           miRNA = mi,
                           bct_mean_fstats = mean(bct$fstats),
                           bmt_mean_fstats = mean(bmt$fstats),
                           bct_mean_pval = mean(bct$pval),
                           bmt_mean_pval = mean(bmt$pval),
                           bct_mean_p.adj = mean(bct$p.adj),
                           bmt_mean_p.adj = mean(bmt$p.adj))
          gene_correlation <- rbind(gene_correlation, dt)
        }
      }
    }
    return(gene_correlation)
  }
  plot.bcmtpColumn <- function(bcmtp, current.col){
    bct <- bcmtp[bcmtp$binding_site_event == "BCT", ]
    bmt <- bcmtp[bcmtp$binding_site_event == "BMT", ]
    par(mar=c(3,5,4,1))
    
    boxplot(bct[[current.col]], bmt[[current.col]],
                       main = paste("Distribution of BCT and BMT correlation ", current.col, sep = ""),
                       ylab = current.col, xlab = "",
                       names = c("bct", "bmt"), las=1)
    cat(paste("Summary of BCT ", current.col, ":\n", sep = ""))
    show(summary(bct[[current.col]]))
    cat(paste("\nSummary of BMT ", current.col, ":\n", sep = ""))
    show(summary(bmt[[current.col]]))
    cat("\n\n")
  }
  
      if(log){
        cat("______________________________________________________________________\n\n")
        cat(paste("\U1F17C\U1F178\U1F181\U1F181\U1F174\U1F183\U1F178", "STATISTICAL ANALYSIS REPORT\n\n\n"))
      }
  gene_correlation <- determin.geneMeans(bcmtp)
  for(current.col in c('coefficient', 'fstats', 'pval', 'p.adj')){
    if(!current.col %in% colnames(bcmtp)){
      next
    }
    cat(paste("Analysis report of correlation ", current.col, "\n", sep = ""))
    cat("-------------------------------------------\n\n")
          
    plot.bcmtpColumn(bcmtp = bcmtp, 
                     current.col = current.col)
    
    #cat("\U1F166\U1F158\U1F15B\U1F152\U1F15E\U1F167\U1F15E\U1F15D-\U1F163\U1F154\U1F162\U1F163\n")
    cat("\U1F146\U1F138\U1F13B\U1F132\U1F13E\U1F147\U1F13E\U1F13D-\U1F143\U1F134\U1F142\U1F143\n")
    test <- wilcox.test(gene_correlation[[paste("bct_mean_", current.col, sep = "")]], 
                        gene_correlation[[paste("bmt_mean_", current.col, sep = "")]],
                        paired = TRUE)
    cat(paste("  V = ", test$statistic, "\n", sep = ""))
    cat(paste("  p-value = ", test$p.value, "\n", sep = ""))
    cat("  alternative hypothesis: true location shift is not equal to 0\n\n")
    
    #cat("\U1F15C\U1F150\U1F15D\U1F15D-\U1F166\U1F157\U1F158\U1F163\U1F15D\U1F154\U1F168-\U1F164-\U1F163\U1F154\U1F162\U1F163\n")
    cat("\U1F13C\U1F130\U1F13D\U1F13D-\U1F146\U1F137\U1F138\U1F143\U1F13D\U1F134\U1F148-\U1F144-\U1F143\U1F134\U1F142\U1F143\n")
    test <- wilcox.test(bcmtp[[current.col]] ~ bcmtp$binding_site_event)
    cat(paste("  W = ", test$statistic, "\n", sep = ""))
    cat(paste("  p-value = ", test$p.value, "\n", sep = ""))
    cat("  alternative hypothesis: true location shift is not equal to 0\n\n")
    
    #cat("\U1F15F\U1F150\U1F158\U1F161\U1F154\U1F153 \U1F163-\U1F163\U1F154\U1F162\U1F163\n")
    cat("\U1F13F\U1F130\U1F138\U1F141\U1F134\U1F133 \U1F143-\U1F143\U1F134\U1F142\U1F143\n")
    test <- t.test(gene_correlation[[paste("bct_mean_", current.col, sep = "")]], 
                   gene_correlation[[paste("bmt_mean_", current.col, sep = "")]],
                   paired = TRUE)
    cat(paste("  t = ", test$statistic, "\n", sep = ""))
    cat(paste("  df = ", test$parameter, "\n", sep = ""))
    cat(paste("  p-value = ", test$p.value, "\n", sep = ""))
    cat("  alternative hypothesis: true location shift is not equal to 0\n")
    cat(paste("  95% confidence interval: ", test$conf.int[1], "\t", test$conf.int[2], "\n", sep = ""))
    cat(paste("  standard error = ", test$stderr, "\n", sep = ""))
    cat("  sample estimates:\n")
    estimate <- as.data.frame(test$estimate)
    for(i in 1:nrow(estimate)){
      cat(paste("    ", rownames(estimate)[i], " = ", estimate[i], "\n", sep = ""))
    }
    cat("\n")
    
    
    #cat("\U1F166\U1F154\U1F15B\U1F152\U1F157 \U1F163\U1F166\U1F15E \U1F162\U1F150\U1F15C\U1F15F\U1F15B\U1F154 \U1F163-\U1F163\U1F154\U1F162\U1F163\n")
    cat("\U1F146\U1F134\U1F13B\U1F132\U1F137 \U1F143\U1F146\U1F13E \U1F142\U1F130\U1F13C\U1F13F\U1F13B\U1F134 \U1F143-\U1F143\U1F134\U1F142\U1F143\n")
    test <- (t.test(bcmtp[[current.col]] ~ bcmtp$binding_site_event))
    cat(paste("  t = ", test$statistic, "\n", sep = ""))
    cat(paste("  df = ", test$parameter, "\n", sep = ""))
    cat(paste("  p-value = ", test$p.value, "\n", sep = ""))
    cat("  alternative hypothesis: true location shift is not equal to 0\n")
    cat(paste("  95% confidence interval: ", test$conf.int[1], "\t", test$conf.int[2], "\n", sep = ""))
    cat(paste("  standard error = ", test$stderr, "\n", sep = ""))
    cat("  sample estimates:\n")
    estimate <- as.data.frame(test$estimate)
    colnames(estimate) <- "V"
    for(i in 1:nrow(estimate)){
      cat(paste("    ", rownames(estimate)[i], " = ", estimate$V[i], "\n", sep = ""))
    }
    cat("\n\n")
  }
      if(log){
        cat("END OF THE STATISTICAL ANALYSIS REPORT\n")
        cat("______________________________________________________________________\n")
      }
}

mirreti.plot.transcriptStatistics <- function(interaction){
  
  determin.transcriptUsage <- function(interactions){
    interactions_corr <- interactions[interactions$correlation != '-', ]
    transcript_bs <- data.table(matrix(nrow = 0, ncol = 5))
    colnames(transcript_bs) <- c('ensembl_gene_id', 'ensembl_transcript_id', 'miRNA', 'bs_numb', 'mean_coefficient')
    for(t in unique(interactions_corr$ensembl_transcript_id)){
      t_sub <- interactions_corr[interactions_corr$ensembl_transcript_id == t, ]
      dt <- data.table(ensembl_gene_id = unique(t_sub$ensembl_gene_id),
                       ensembl_transcript_id = t)
    }
  }
}

  

# #' DESCRIBE METHODE
# #'
# #' @param tpm_expr
# #' @param mir_expr
# #' @param sample_annotation
# #' @param interactions
# #' @param cluster.size
# #' @param log
# #' @return interactions   #, drimseq
mirreti.twoConditions <- function(tpm_expr, mir_expr, sample_annotation, interactions, cluster.size, log = FALSE){
  stopifnot(uniqueN(sample_annotation$condition) == 2)
  
  #dubmethods
  extractSamples.forCondition <- function(tpm_expr, mir_expr, sample_annotation, condition){
    sample_annotation <- sample_annotation[sample_annotation$condition == condition, ]
    tpm_expr <- tpm_expr[ , colnames(tpm_expr) %in% sample_annotation$sample]
    mir_expr <- mir_expr[ , colnames(mir_expr) %in% sample_annotation$sample]
    return(list(tpm_expr, mir_expr))
  }
  
  
  # Step 1: run SPONGE interaction filter for each condition in sample_annotation
  for(c in unique(sample_annotation$condition)){
    expr_data <- extractSamples.forCondition(tpm_expr = tpm_expr,
                                             mir_expr = mir_expr,
                                             sample_annotation = sample_annotation,
                                             condition = c)
    interactions <- sponge.correlateData(tpm_expr = expr_data[[1]],
                                         mir_expr = expr_data[[2]],
                                         interactions = interactions,
                                         condition = c,
                                         cluster.size = cluster.size,
                                         log = log)
  }
 
  
  # Step 2: determin DTU with DRIMSeq
  #drimseq <- drimseq.determin.dtu(tpm_expr = tpm_expr,
  #                                sample_annotation = sample_annotation,
  #                                log = TRUE)
  
  #return(list(interactions, drimseq))
  return(interactions)
}




#------------------------------------------------------------------------------------------------
#                                     SPONGE related methods
#                              - used directly in MIRRETI methods -
#------------------------------------------------------------------------------------------------

# #' DESCRIBE METHODE
# #'
# #' @param tpm_expr
# #' @param mir_expr
# #' @param interactions
# #' @param condition
# #' @param f.test
# #' @param cluster.size
# #' @param log
# #' @return interactions
sponge.correlateData <- function(tpm_expr, mir_expr, interactions, condition = NULL, f.test = FALSE, cluster.size, log = FALSE){
      start.time <- 0
      end.time <- 0
  
      
  # submethods
  unlist.candidates <- function(candidates){
    candidates <- do.call(rbind, candidates)
    candidates$transcript <- lapply(strsplit(rownames(candidates), "\\."), '[', 1)
    rownames(candidates) <- NULL
    return(candidates)
  }
  summarize.candidates <- function(candidates, colname){
        column <- dplyr::select(candidates, all_of(colname))
        summary <- summary(column)
        cat(paste("    summary of correlation ", colname, "\n", sep = ""))
        cat(paste("     - Min.: ", summary[1], "\n", sep = ""))
        cat(paste("     - 1st Qu.: ", summary[2], "\n", sep = ""))
        cat(paste("     - Median: ", summary[3], "\n", sep = ""))
        cat(paste("     - Mean: ", summary[4], "\n", sep = ""))
        cat(paste("     - 3rd Qu.: ", summary[5], "\n", sep = ""))
        cat(paste("     - Max.: ", summary[6], "\n", sep = ""))
        boxplot(column, ylab = colname, xlab = "interactions", main = paste("Distribution of the negativly correlated\ninteraction ", colname, sep = ""))
  }
  integrate.candidates.intoInteractions <- function(interactions, candidates){
    interactions_mrnamirna_pairs <- as.character(unique(paste(interactions$miRNA, interactions$ensembl_transcript_id, sep = ",")))
    interactions_correlated <- as.character(unique(paste(candidates$miRNA, candidates$ensembl_transcript_id, sep = ",")))
    interactions_notCorrelated <- interactions_mrnamirna_pairs[!interactions_mrnamirna_pairs %in% interactions_correlated]
      rm(interactions_mrnamirna_pairs, interactions_correlated)
      gc()
    if(ncol(candidates) == 4){
      interactions_notCorrelated <- data.table(miRNA = lapply(strsplit(interactions_notCorrelated, ','), '[', 1),
                                               ensembl_transcript_id = lapply(strsplit(interactions_notCorrelated, ','), '[', 2),
                                               coefficient = rep(0, length(interactions_notCorrelated)),
                                               correlation = rep("-", length(interactions_notCorrelated)))
    }
    else if(ncol(candidates) == 6){
      colnames(candidates)
      interactions_notCorrelated <- data.table(miRNA = lapply(strsplit(interactions_notCorrelated, ','), '[', 1),
                                               ensembl_transcript_id = lapply(strsplit(interactions_notCorrelated, ','), '[', 2),
                                               correlation = rep("-", length(interactions_notCorrelated)),
                                               fstats = rep(0, length(interactions_notCorrelated)),
                                               pval = rep(0, length(interactions_notCorrelated)),
                                               p.adj = rep(0, length(interactions_notCorrelated)))
    }
    else{
      cat("candidates colnames:\n")
      cat(colnames(candidates))
      cat("\n")
      stopifnot(FALSE)
    }
    correlation_table <- rbind(candidates, interactions_notCorrelated)
      rm(interactions_notCorrelated)
      gc()
    correlation_table$miRNA <- as.character(correlation_table$miRNA)
    correlation_table$ensembl_transcript_id <- as.character(correlation_table$ensembl_transcript_id)
    interactions <- merge(interactions, correlation_table, by = c("miRNA", "ensembl_transcript_id"))
    return(interactions)
  }
  
  
  #Step 1: process data for SPONGE interaction filter
      if(log){
        cat("______________________________________________________________________\n\n")
        cat(paste("\U1F182\U1F17F\U1F17E\U1F17D\U1F176\U1F174", "DATA CORRELATION REPORT\n"))
        cat("\n  \U1F664 \U1F158\U1F15D\U1F15F\U1F164\U1F163 \U1F161\U1F154\U1F15F\U1F15E\U1F161\U1F163 \U1F666\n")
        cat("  INTERACTIONS\n")
        cat(paste("    Total interactions: ", nrow(interactions), "\n", sep = ""))
        cat(paste("    transcripts: ", uniqueN(interactions$ensembl_transcript_id), "\n", sep = ""))
        cat(paste("    genes: ", uniqueN(interactions$ensembl_gene_id), "\n", sep = ""))
        cat(paste("    microRNAs: ", uniqueN(interactions$miRNA), "\n", sep = ""))
        cat("  TRANSCRIPT EXPRESSION DATA\n")
        cat(paste("    transcripts: ", uniqueN(rownames(tpm_expr)), "\n", sep = ""))
        cat(paste("    number of samples: ", uniqueN(colnames(tpm_expr)), "\n", sep = ""))
        cat("  MICRORNA EXPRESSION DATA\n")
        cat(paste("    microRNAs: ", uniqueN(rownames(mir_expr)), "\n", sep = ""))
        cat(paste("    number of samples: ", uniqueN(colnames(mir_expr)), "\n", sep = ""))
        cat(paste("  CONDITION: ", condition, "\n", sep = ""))
        cat(paste("  F-TEST: ", f.test, "\n", sep = ""))
        cat(paste("  CLUSTER SIZE: ", cluster.size, "\n\n\n", sep = ""))
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart preparing expression and interactions data for SPONGE interaction filter...\n", sep = ""))
      }
  mir_expr <- t(mir_expr)
  tpm_expr <- t(tpm_expr)
  interactions_matrix <- as.matrix(as.data.frame.matrix(table(interactions[ , c("ensembl_transcript_id", "miRNA")])))
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\t\U2705 \U1F1F8\U1F1FA\U1F1E8\U1F1E8\U1F1EA\U1F1F8\U1F1F8\U1F1EB\U1F1FA\U1F1F1\n", sep = ""))
        cat(paste("\U231B Time required:\t", round(difftime(end.time, start.time, units = "secs"), digits = 3), " secs\n\n", sep = ""))
      } 
  
  
  # Step 2: run SPONGE interaction filter
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart running SPONGE interaction filter...\n", sep = ""))
      }
    cl <- makePSOCKcluster(cluster.size)
    registerDoParallel(cl)
  candidates <- sponge_gene_miRNA_interaction_filter(gene_expr = tpm_expr,
                                                     mir_expr = mir_expr,
                                                     mir_predicted_targets = interactions_matrix,
                                                     coefficient.threshold = -0.05,
                                                     coefficient.direction = "<",
                                                     F.test = f.test,
                                                     parallel.chunks = cluster.size)
    stopCluster(cl)
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\t\U2705 \U1F1F8\U1F1FA\U1F1E8\U1F1E8\U1F1EA\U1F1F8\U1F1F8\U1F1EB\U1F1FA\U1F1F1\n", sep = ""))
        cat(paste("\U231B Time required:\t", round(difftime(end.time, start.time, units = "secs"), digits = 3), " secs\n\n", sep = ""))
      }
  
  
  # Step 3: unlist SPONGE candidates
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart integrating correlation coefficient and direction of SPONGE candidates into interactions data...\n", sep = ""))
      }
  candidates <- unlist.candidates(candidates = candidates)
  candidates$correlation <- rep("<", nrow(candidates))
  if(!f.test){
    colnames(candidates) <- c('miRNA', 'coefficient', 'ensembl_transcript_id', 'correlation')
  }
  else{
    colnames(candidates) <- c('miRNA', 'fstats', 'pval', 'p.adj', 'ensembl_transcript_id', 'correlation')
  }
  
  
  # Step 4: integrate coefficients and correlation information into interaction data
  interactions <- integrate.candidates.intoInteractions(interactions, candidates)
  candidates <- interactions[interactions$correlation != '-', ]
  if(!is.null(condition) & 'coefficient' %in% colnames(interactions)){
    colnames(interactions)[which(colnames(interactions) == "correlation")] <- paste('correlation', condition, sep = "_")
    colnames(interactions)[which(colnames(interactions) == "coefficient")] <- paste('coefficient', condition, sep = "_")
  }
  else if(!is.null(condition) & 'fstats' %in% colnames(interactions) & 'pval' %in% colnames(interactions) & 'p.adj' %in% colnames(interactions)){
    colnames(interactions)[which(colnames(interactions) == "correlation")] <- paste('correlation', condition, sep = "_")
    colnames(interactions)[which(colnames(interactions) == "fstats")] <- paste('fstats', condition, sep = "_")
    colnames(interactions)[which(colnames(interactions) == "pval")] <- paste('pval', condition, sep = "_")
    colnames(interactions)[which(colnames(interactions) == "p.adj")] <- paste('p.adj', condition, sep = "_")
  }
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\t\U2705 \U1F1F8\U1F1FA\U1F1E8\U1F1E8\U1F1EA\U1F1F8\U1F1F8\U1F1EB\U1F1FA\U1F1F1\n", sep = ""))
        cat(paste("\U231B Time required:\t", round(difftime(end.time, start.time, units = "secs"), digits = 3), " secs\n\n", sep = ""))
        cat("\n  \U1F664 \U1F15E\U1F164\U1F163\U1F15F\U1F164\U1F163 \U1F161\U1F154\U1F15F\U1F15E\U1F161\U1F163 \U1F666\n")
        cat("  INTERACTIONS\n")
        cat(paste("    Total interactions: ", nrow(interactions), "\n", sep = ""))
        cat(paste("    correlated interactions: ", nrow(candidates), "\n", sep = ""))
        cat(paste("    correlated transcripts: ", uniqueN(candidates$ensembl_transcript_id), "\n", sep = ""))
        cat(paste("    correlated genes: ", uniqueN(candidates$ensembl_gene_id), "\n", sep = ""))
        cat(paste("    correlated microRNAs: ", uniqueN(candidates$miRNA), "\n", sep = ""))
        for(colname in c('coefficient', 'fstats', 'pval', 'p.adj')){
          if(!colname %in% colnames(candidates)){
            next
          }
          summarize.candidates(candidates = candidates,
                               colname = colname)
        }
        cat("\n\nEND OF THE SPONGE DATA CORRELATION REPORT\n")
        cat("______________________________________________________________________\n")
      } 
  
  
  return(interactions)
}



#------------------------------------------------------------------------------------------------
#                                     DRIMSeq related methods
#                              - used directly in MIRRETI methods -
#------------------------------------------------------------------------------------------------

# #' DESCRIBE METHODE
# #'
# #' @param tpm_expr
# #' @param sample_annotation
# #' @param ensembl_mart
# #' @param log
# #' @return drimseq
drimseq.determin.dtu <- function(tpm_expr, sample_annotation, ensembl_mart = NULL, log = FALSE){
      start.time <- 0
      end.time <- 0
  
  
  # Step 1: build gene model from Ensembl
      if(log){
        cat("\n\nMIRRETI DRIMSeq report\n")
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart building gene model on Ensembl transcripts...\n", sep = ""))
      }
  if(is.null(ensembl_mart)){
    ensembl_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  }
  biomart_annot <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"),
                         filters = "ensembl_transcript_id",
                         values = rownames(tpm_expr),
                         mart = ensembl_mart)
  biomart_annot <- dplyr::select(biomart_annot, ensembl_gene_id, ensembl_transcript_id)
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
  
  
  # Step 2: prepare samples and counts data for DEA
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart preparing samples and counts data for differential expression...\n", sep = ""))
      }
  tpm_expr <- tpm_expr[rownames(tpm_expr) %in% biomart_annot_multitranscript$feature_id, ]
  tpm_expr <- 2^tpm_expr
  tpm_expr$feature_id <- rownames(tpm_expr)
  tpm_expr <- cbind(biomart_annot_multitranscript, tpm_expr, by = "feature_id")
  rownames(tpm_expr) <- NULL
  sample_annotation <- sample_annotation[sample_annotation$sample %in% colnames(tpm_expr), ]
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\tPreparing data for DRIMSeq was successfull\n", sep = ""))
        cat(paste("Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
      }
  
  
  # Step 3: DRIMSeq
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart running DRIMSeq...\n", sep = ""))
      }
  drimseq <- dmDSdata(counts=tpm_expr,
                      samples=sample_annotation)
  #table(table(counts(d)$gene_id))
  design_full <- model.matrix(~condition, data=DRIMSeq::samples(drimseq))
  #colnames(design_full)
  
    cl <- makePSOCKcluster(cluster.size)
    registerDoParallel(cl)
  system.time({
    drimseq <- dmPrecision(drimseq, design=design_full)
    drimseq <- dmFit(drimseq, design=design_full)
    drimseq <- dmTest(drimseq, coef="condition")
  })
    stopCluster(cl)
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\tRunning DRIMSeq was successfull\n", sep = ""))
        cat(paste("Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
        cat("End of the MIRRETI DRIMSeq report")
      }
  
  return(drimseq)
}

