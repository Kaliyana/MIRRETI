#load("~/bachelor/data/core_data/mirreti_coredata.RDS")
#set.seed(01101990)
library(data.table)                                      # uniqueN(), unique()
library(dplyr)                                           # %>%, select()
library(crayon)                                          # crayon::bold()
library(biomaRt)                                         # useMart(), getBM(), listMarts()
library(doParallel)                                      # makePSOCKcluster(), registerDoParallel(), stopCluster()
library(SPONGE)                                          # sponge_gene_miRNA_interaction_filter()




# #' DESCRIBE METHODE
# #'
# #' @param tpm_expr
# #' @param mir_expr
# #' @param sample_annotation
# #' @param primary.diseases
# #' @param conditions
# #' @param log
# #' @return tpm_expr, mir_expr, sample_annotation
tcgaSamples.subset <- function(tpm_expr, mir_expr, sample_annotation, primary.diseases, conditions = "1", log = FALSE){
      start.time <- 0
      end.time <- 0
  
  # Step 1: filter samples 
      if(log){
        cat(paste("\n\n\U1F183\U1F172\U1F176\U1F170", crayon::bold("expression data subset report\n")))
        cat("\n  \U1F664 \U1F158\U1F15D\U1F15F\U1F164\U1F163 \U1F161\U1F154\U1F15F\U1F15E\U1F161\U1F163 \U1F666\n")
        cat(crayon::bold("  samples\n"))
        cat(paste("\tTotal samples: ", uniqueN(sample_annotation$sample), "\n", sep = ""))
        cat(paste("\tNumber of cancer types: ", uniqueN(sample_annotation$X_primary_disease), "\n", sep = ""))
        cat(crayon::bold("  transcript expression data\n"))
        cat(paste("\tTranscripts: ", uniqueN(rownames(tpm_expr)), "\n", sep = ""))
        cat(paste("\tNumbeer of samples: ", uniqueN(colnames(tpm_expr)), "\n", sep = ""))
        cat(crayon::bold("  microRNA expression data\n"))
        cat(paste("\tmicroRNAs: ", uniqueN(rownames(mir_expr)), "\n", sep = ""))
        cat(paste("\tNumber of samples: ", uniqueN(colnames(mir_expr)), "\n", sep = ""))
        cat(paste("\tcancer types: ", paste(primary.diseases, collapse = "|"), "\n", sep = ""))
        cat(paste("\tconditions: ", paste(conditions, collapse = "|"), "\n\n", sep = ""))
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart filtering TCGA samples for cancer types and conditions...\n", sep = ""))
      }
  sample_annotation <- sample_annotation %>% filter(X_primary_disease %in% primary.diseases &
                                                      sample_type_id %in% conditions)
  sample_annotation$condition <- paste(gsub(" ", "", sample_annotation$X_primary_disease), sample_annotation$sample_type_id, sep = "_")
  sample_annotation <- dplyr::select(sample_annotation, sample, condition)
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, crayon::bold("\t\U1F5F8 \U1F1F8\U1F1FA\U1F1E8\U1F1E8\U1F1EA\U1F1F8\U1F1F8\U1F1EB\U1F1FA\U1F1F1\n"), sep = ""))
        cat(paste("\tTotal samples: ", uniqueN(sample_annotation$sample), "\n", sep = ""))
        cat("\t!!!!!!!!Hier könnten Ihre Frequenzen stehen!!!!!!!!\n")
        cat(paste("Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
      } 
  
  
  # Step 2: subset TPM expression data and filter transcripts and samples for NAs
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart subsetting transcript expression data for samples and filter for NAs...\n", sep = ""))
      }
  tpm_expr <- tpm_expr[ , colnames(tpm_expr) %in% sample_annotation$sample]
  tpm_expr <- exprdata.filter.na(tpm_expr)
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, crayon::bold("\t\U1F5F8 \U1F1F8\U1F1FA\U1F1E8\U1F1E8\U1F1EA\U1F1F8\U1F1F8\U1F1EB\U1F1FA\U1F1F1\n"), sep = ""))
        cat(paste("\tTranscripts: ", uniqueN(rownames(tpm_expr)), "\n", sep = ""))
        cat(paste("\tRemaining samples: ", uniqueN(colnames(tpm_expr)), "\n", sep = ""))
        cat(paste("Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
      }
  
  
  # Step 3: subset miRNA expression data and filter microRNAs and samples for NAs
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart subsetting microRNA expression data for samples and filter for NAs...\n", sep = ""))
      }
  mir_expr <- mir_expr[ , colnames(mir_expr) %in% sample_annotation$sample]
  mir_expr <- exprdata.filter.na(mir_expr)
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, crayon::bold("\t\U1F5F8 \U1F1F8\U1F1FA\U1F1E8\U1F1E8\U1F1EA\U1F1F8\U1F1F8\U1F1EB\U1F1FA\U1F1F1\n"), sep = ""))
        cat(paste("\tmicroRNAs: ", uniqueN(rownames(mir_expr)), "\n", sep = ""))
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
        cat(paste(end.time, crayon::bold("\t\U1F5F8 \U1F1F8\U1F1FA\U1F1E8\U1F1E8\U1F1EA\U1F1F8\U1F1F8\U1F1EB\U1F1FA\U1F1F1\n"), sep = ""))
        cat(paste("Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
        cat("\n  \U1F664 \U1F15E\U1F164\U1F163\U1F15F\U1F164\U1F163 \U1F161\U1F154\U1F15F\U1F15E\U1F161\U1F163 \U1F666\n")
        cat(crayon::bold("  samples\n"))
        cat(paste("\tTotal samples: ", uniqueN(sample_annotation$sample), "\n", sep = ""))
        cat("\t!!!!!!!!Hier könnten Ihre Frequenzen stehen!!!!!!!!\n")
        cat(crayon::bold("  transcript expression data\n"))
        cat(paste("\tTranscripts: ", uniqueN(rownames(tpm_expr)), "\n", sep = ""))
        cat(paste("\tRemaining samples: ", uniqueN(colnames(tpm_expr)), "\n", sep = ""))
        cat(crayon::bold("  microRNA expression data\n"))
        cat(paste("\tmicroRNAs: ", uniqueN(rownames(mir_expr)), "\n", sep = ""))
        cat(paste("\tRemaining samples: ", uniqueN(colnames(mir_expr)), "\n\n", sep = ""))
        cat(crayon::bold("End of the TCGA expression data subset report\n\n"))
      }
  
  return(list(tpm_expr, mir_expr, sample_annotation))
}

exprdata.filter.na <- function(expr_data, max.na.percent = 0.2){
  na_content <- rowSums(expr_data == min(expr_data))
  expr_data <- expr_data[na_content/ncol(expr_data) <= max.na.percent, ]
  na_content <- colSums(expr_data == min(expr_data))
  expr_data <- expr_data[ , na_content/nrow(expr_data) <= max.na.percent]
  return(expr_data)
}



# #' Search RefSeq mRNA IDs of your miRWalk interactions data in Ensembl database to obtain matching Ensembl transcript IDs.
# #'
# #' @param interactions: Data Frame of miRNA - mRNA interactions. The column name of the mRNAs needs to be "mRNA".
# #' @param ensembl_mart: Mart object compatible with the biomaRt query. Gives the possibility to search older versions of the Ensembl database. The default is the current Ensembl database version.
# #' @param id.type: Search Ensembl database for "validated", "predicted" or "both" RefSeq mRNA IDs. The default will search for both ID types.
# #' @return interactions: Modified interactions data with matching Ensembl transcript IDs
mirwalkInteractions.idconversion <- function(interactions, ensembl_mart = NULL, id.type = 'both'){
  
  # Step 1: Prepare reqired data
  refseq_ids <- unique(interactions$mRNA)
  if(is.null(ensembl_mart))
    ensembl_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  
  # Step 2: Retrieve Ensemble IDs via biomaRt
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
  }
  else {
    cat("\nIncorrect value for parameter id_type\nPlease choose between 'validated', 'predicted' and 'both'\nRefSeq to Ensemble IDs conversion stopped...\n")
    stopifnot(FALSE)
  }
  
  stopifnot(!is.null(mappings))
  
  
  # Step 3: Remove entries where no Ensemble IDs were found
  # THIS MIRRETI VERSION USES ONLY DISTINCT MAPPINGS  
  ambiguous_mrnas <- mappings[duplicated(mappings$mRNA), ]
  ambiguous_transcripts <- mappings[duplicated(mappings$ensembl_transcript_id), ]
  distinct_mappings <- mappings[!(mappings$mRNA %in% ambiguous_mrnas$mRNA) &
                                  !(mappings$ensembl_transcript_id %in% ambiguous_transcripts$ensembl_transcript_id), ]
  
  distinct_interactions <- interactions[interactions$mRNA %in% distinct_mappings$mRNA, ]
  interactions <- merge(distinct_interactions, distinct_mappings, by = c("mRNA"), allow.cartesian = TRUE)
  interactions <- dplyr::select(interactions, all_of(additional.attributes), mRNA, miRNA, binding_site, binding_probability, genomic_region)
  
  return(interactions)
    rm(refseq_ids, mappings, refseq_ensembl, predrefseq_ensembl, ambiguous_mrnas, ambiguous_transcripts, distinct_interactions, distinct_mappings, id.type)
}



# #' DESCRIBE METHODE
# #'
# #' @param tpm_expr
# #' @param mir_expr
# #' @param interactions
# #' @param var.threshold
# #' @return tpm_expr, mir_expr, interactions
mirretiData.filter.actualTargets <- function(tpm_expr, mir_expr, interactions, var.threshold = 0.2){
  
  # Step 1: filter out var.threashold least variant columns (molecules) in expression data
  rownames(tpm_expr) <- lapply(strsplit(rownames(tpm_expr), "\\."), '[', 1)
  tpm_expr <- exprdata.filter.variance(expr_data = tpm_expr,
                                       var.threshold = var.threshold)
  mir_expr <- exprdata.filter.variance(expr_data = mir_expr,
                                       var.threshold = var.threshold)
  
  # Step 2: match actual overlapping microRNAs and transcript in expression and interactions data
  tpm_expr <- tpm_expr[rownames(tpm_expr) %in% interactions$ensembl_transcript_id, colnames(tpm_expr) %in% colnames(mir_expr)]
  mir_expr <- mir_expr[rownames(mir_expr) %in% interactions$miRNA , colnames(mir_expr) %in% colnames(tpm_expr)]
  interactions <- interactions[interactions$miRNA %in% rownames(mir_expr) & interactions$ensembl_transcript_id %in% rownames(tpm_expr), ]

  return(list(tpm_expr, mir_expr, interactions))
}

exprdata.filter.variance <- function(expr_data, var.threshold = 0.2){
  expr_data$variance <- apply(data.frame(expr_data), 1, var)
  expr_data <- expr_data[order(expr_data$variance, decreasing = T), ]
  threshold <- floor((1-var.threshold) * nrow(expr_data))
  expr_data <- expr_data[1:threshold, colnames(expr_data) != 'variance']
  return(expr_data)
}



# #' DESCRIBE METHODE
# #'
# #' @param tpm_expr
# #' @param mir_expr
# #' @param interactions
# #' @param cluster.size
# #' @param log
# #' @return interactions
mirreti.simplerun <- function(tpm_expr, mir_expr, interactions, cluster.size, log = FALSE){
      start.time <- 0
      end.time <- 0
  
  
  #Step 1: process data for SPONGE interaction filter
      if(log){
        cat(paste("\n\n\U1F17C\U1F178\U1F181\U1F181\U1F174\U1F183\U1F178", crayon::bold("simple run report\n")))
        cat("\n  \U1F664 \U1F158\U1F15D\U1F15F\U1F164\U1F163 \U1F161\U1F154\U1F15F\U1F15E\U1F161\U1F163 \U1F666\n")
        cat(crayon::bold("  interactions\n"))
        cat(paste("\tTotal interactions: ", nrow(interactions), "\n", sep = ""))
        cat(paste("\ttranscripts: ", uniqueN(interactions$ensembl_transcript_id), "\n", sep = ""))
        cat(paste("\tgenes: ", uniqueN(interactions$ensembl_gene_id), "\n", sep = ""))
        cat(paste("\tmicroRNAs: ", uniqueN(interactions$miRNA), "\n", sep = ""))
        cat(crayon::bold("  transcript expression data\n"))
        cat(paste("\tTranscripts: ", uniqueN(rownames(tpm_expr)), "\n", sep = ""))
        cat(paste("\tNumbeer of samples: ", uniqueN(colnames(tpm_expr)), "\n", sep = ""))
        cat(crayon::bold("  microRNA expression data\n"))
        cat(paste("\tmicroRNAs: ", uniqueN(rownames(mir_expr)), "\n", sep = ""))
        cat(paste("\tNumber of samples: ", uniqueN(colnames(mir_expr)), "\n", sep = ""))
        cat(paste("\tcluster size: ", cluster.size, "\n\n", sep = ""))
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart preparing expression and interactions data for SPONGE interaction filter...\n", sep = ""))
      }
  mir_expr <- t(mir_expr)
  tpm_expr <- t(tpm_expr)
  interactions_matrix <- as.matrix(as.data.frame.matrix(table(interactions[ , c("ensembl_transcript_id", "miRNA")])))
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, crayon::bold("\t\U1F5F8 \U1F1F8\U1F1FA\U1F1E8\U1F1E8\U1F1EA\U1F1F8\U1F1F8\U1F1EB\U1F1FA\U1F1F1\n"), sep = ""))
        cat(paste("Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
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
                                                     parallel.chunks = cluster.size)
    stopCluster(cl)
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, crayon::bold("\t\U1F5F8 \U1F1F8\U1F1FA\U1F1E8\U1F1E8\U1F1EA\U1F1F8\U1F1F8\U1F1EB\U1F1FA\U1F1F1\n"), sep = ""))
        cat(paste("Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
      }

  
  # Step 3: integrate coefficients and correlation information into interaction data
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart integrating correlation coefficient and direction of SPONGE candidates into interactions data...\n", sep = ""))
      }
  candidates <- sponge.candidates.unlist(candidates)
  candidates$correlation <- rep("<", nrow(candidates))
  colnames(candidates) <- c('miRNA', 'coefficient', 'ensembl_transcript_id', 'correlation')
  interactions <- sponge.integrate.correlatedInteractions(interactions, candidates)
  candidates <- interactions[interactions$correlation == "<", ]
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, crayon::bold("\t\U1F5F8 \U1F1F8\U1F1FA\U1F1E8\U1F1E8\U1F1EA\U1F1F8\U1F1F8\U1F1EB\U1F1FA\U1F1F1\n"), sep = ""))
        cat(paste("Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
        cat("\n  \U1F664 \U1F15E\U1F164\U1F163\U1F15F\U1F164\U1F163 \U1F161\U1F154\U1F15F\U1F15E\U1F161\U1F163 \U1F666\n")
    cat(crayon::bold("  interactions\n"))
    cat(paste("\tTotal interactions: ", nrow(interactions), "\n", sep = ""))
    cat(paste("\tcorrelated interactions: ", nrow(candidates), "\n", sep = ""))
    cat(paste("\tcorrelated transcripts: ", uniqueN(candidates$ensembl_transcript_id), "\n", sep = ""))
    cat(paste("\tcorrelated genes: ", uniqueN(candidates$ensembl_gene_id), "\n", sep = ""))
    cat(paste("\tcorrelated microRNAs: ", uniqueN(candidates$miRNA), "\n", sep = ""))
    cat(paste("\tsummary of correlation coefficients: ", summary(candidates$coefficient), "\n", sep = ""))
    boxplot(correlation_table$coefficient, ylab = "coefficient", xlab = "interactions", main = "Distribution of the negativly correlated\ninteraction coefficients")
    cat(crayon::bold("\n\nEnd of the SPONGE simple run report\n\n"))
  }
  
  return(interactions)
}

sponge.candidates.unlist <- function(candidates){
  candidates <- do.call(rbind, candidates)
  candidates$transcript <- lapply(strsplit(rownames(candidates), "\\."), '[', 1)
  rownames(candidates) <- NULL
  return(candidates)
}

sponge.integrate.correlatedInteractions <- function(interactions, candidates){
  interactions_mrnamirna_pairs <- as.character(unique(paste(interactions$miRNA, interactions$ensembl_transcript_id, sep = ",")))
  interactions_correlated <- as.character(unique(paste(candidates$miRNA, candidates$ensembl_transcript_id, sep = ",")))
  interactions_notCorrelated <- interactions_mrnamirna_pairs[!interactions_mrnamirna_pairs %in% interactions_correlated]
      rm(interactions_mrnamirna_pairs, interactions_correlated)
      gc()
  interactions_notCorrelated <- data.table(miRNA = lapply(strsplit(interactions_notCorrelated, ','), '[', 1),
                                           ensembl_transcript_id = lapply(strsplit(interactions_notCorrelated, ','), '[', 2),
                                           coefficient = rep(0, length(interactions_notCorrelated)),
                                           correlation = rep("-", length(interactions_notCorrelated)))
  correlation_table <- rbind(candidates, interactions_notCorrelated)
  rm(interactions_notCorrelated)
  correlation_table$miRNA <- as.character(correlation_table$miRNA)
  correlation_table$ensembl_transcript_id <- as.character(correlation_table$ensembl_transcript_id)
  interactions <- merge(interactions, correlation_table, by = c("miRNA", "ensembl_transcript_id"))
  return(interactions)
}
