library(data.table)
library(graphics)
library(tidyverse)
library(dplyr)
library(biomaRt)
library(SPONGE)
library(doParallel)
#set.seed(01101990)
#ensembl_mart_99 <- useMart("ENSEMBL_MART_ENSEMBL", host = "http://jan2020.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")
#ensembl_mart_92 <- useMart("ENSEMBL_MART_ENSEMBL", host = "http://apr2018.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")


#-------------------------------------------------------------------------------------
#                                     MIRRETI method
#-------------------------------------------------------------------------------------


mirreti <- function(mir.filepath, tpm.filepath, interactions.filepath, sampleannot.filepath, primary.disease, sample.type.id, cluster.size = 40, top.number = NULL, ensembl_mart = NULL){
  
      cat("\n\t\t\t\t*** Welcome to MIRRETI ***\n\n")
  
  
  # Step 1: Reading in the data
      start.time <- Sys.time()
      cat("______________________________________________________________________\n")
      cat(paste(start.time, "\t\tStart reading the data...\n", sep = ""))
  data <- read.mirreti.data(mir.filepath = mir.filepath,
                            tpm.filepath = tpm.filepath,
                            interactions.filepath = interactions.filepath,
                            sampleannot.filepath = sampleannot.filepath)
  mir_expr <- data[[1]]
  tpm_expr <- data[[2]]
  interactions <- data[[3]]
  sample_annotation <- data[[4]]
      end.time <- Sys.time()
      diff.time <- difftime(end.time, start.time, units = "secs")
      cat(paste(end.time, "\t\tFinished reading the data\n=> Total procession time:\t", diff.time, " SECS\n\n", sep = ""))
  
  
  # Step 2: Preprocessing the data
      start.time <- Sys.time()
      cat("______________________________________________________________________\n")
      cat(paste(start.time, "\t\tStart preprocessing the data...\n", sep = ""))
  data <- preproc(mir_expr = mir_expr,
                  tpm_expr = tpm_expr,
                  interactions = interactions,
                  sample_annotation = sample_annotation,
                  primary.disease = primary.disease,
                  sample.type.id = sample.type.id,
                  ensembl_mart = ensembl_mart, 
                  top.number = NULL)
  mir_expr <- data[[1]]
  tpm_expr <- data[[2]]
  interactions <- data[[3]]
  interactions_matrix <- data[[4]]
      end.time <- Sys.time()
      diff.time <- difftime(end.time, start.time, units = "secs")
      cat(paste(end.time, "\t\tFinished preprocessing the data\n=> Total procession time:\t", diff.time, " SECS\n\n", sep = ""))
  
  # Overview on data dimensions
      cat("Report on data dimensions handed to SPONGE interaction filter\n")
      cat(paste("\t-> miRNA: ", nrow(mir_expr), " samples & ", ncol(mir_expr), " columns\n"))
      cat(paste("\t-> transcripts: ", nrow(tpm_expr), " rows & ", ncol(tpm_expr), " columns\n"))
      cat(paste("\t-> interactions: ", nrow(interactions_matrix), " rows & ", ncol(interactions_matrix), " columns\n\n"))
  
    
  # Step 3: Run SPONGE interaction filter and combine with annotation
  # Notifications in sponge.filter()
  candidates <- sponge.filter(tpm_expr = tpm_expr,
                              mir_expr = mir_expr,
                              interactions_matrix, cluster.size)

  
  #Step 4: Identify exons where candidate binding sites are located
      start.time <- Sys.time()
      cat("______________________________________________________________________\n")
      cat(paste(start.time, "\t\tStart enriching the exon model...\n", sep = ""))
  candidate_events <- exon.model.enrich(interactions = interactions,
                                        candidates = candidates,
                                        ensembl_mart = ensembl_mart_99)
  candidate_events <- candidate_events %>% filter(ensembl_transcript_id_1 %in% colnames(tpm_expr) & transcript_without_exon_and_correlation %in% colnames(tpm_expr)) 
  tpm_expr <- tpm_expr[ , colnames(tpm_expr) %in% candidate_events$ensembl_transcript_id_1 |
                          colnames(tpm_expr) %in% candidate_events$transcript_without_exon_and_correlation]
      end.time <- Sys.time()
      diff.time <- difftime(end.time, start.time, units = "secs")
      cat(paste(end.time, "\t\tFinished enriching the exon model\n=> Total procession time:\t", diff.time, " SECS\n\n", sep = ""))
  
  # Output report
      temp <- unique(c(as.vector(candidate_events$ensembl_transcript_id_1), as.vector(candidate_events$transcript_without_exon_and_correlation)))
      cat("Report on MIRRETI output\n")
      cat(paste("\t", length(temp), " unique transcripts in candidate_events (ensembl_transcript_id_1 & transcript_without_Exon_and_correlation columns)"))
      cat(paste("\t", ncol(tpm_expr), " unique transcripts in tpm_expr"))
      cat(paste("\t", length(intersect(temp, colnames(tpm_expr))), " common transcripts between candidate_events & tpm_expr"))
      
      
      cat("\t*** END OF THE MIRRETI RUNTIME ANALYSIS ***")
  
  return(list(candidate_events, tpm_expr))
}



#-------------------------------------------------------------------------------------
#                                 Basic read methods
#-------------------------------------------------------------------------------------


read.expression.file <- function(expression.filepath){
  expr_data <- data.frame(fread(expression.filepath, header = T, sep = "\t"), row.names = 1, check.names = F)
  rownames(expr_data) <- lapply(strsplit(rownames(expr_data), "\\."), '[', 1)
  return(expr_data)
}

read.interaction.file <- function(interaction.filepath){
  interactions <- fread(interaction.filepath, header = T, sep = "\t")
  return(interactions)
}

read.interaction.files <- function(utr5.filepath, cds.filepath, utr3.filepath){
  utr5 <- read.interaction.file(utr5.filepath)
  utr5$genomic_region <- rep("5UTR", nrow(utr5))
  cds <- read.interaction.file(cds.filepath)
  cds$genomic_region <- rep("CDS", nrow(cds))
  utr3 <- read.interaction.file(utr3.filepath)
  utr3$genomic_region <- rep("3UTR", nrow(utr3))
  interactions <- rbind(utr5, cds, utr3)
  return(interactions)
}

read.sample.annotation.file <- function(sampleannotation.filepath){
  sample_annotation <- read.csv(sampleannotation.filepath, header = T, sep = "\t")
  return(sample_annotation)
}

# MASTER METHOD
read.mirreti.data <- function(mir.filepath, tpm.filepath, interactions.filepath, sampleannot.filepath){
  mir_expr <- read.expression.file(mir.filepath)
  tpm_expr <- read.expression.file(tpm.filepath)
  interactions <- NULL
  if(typeof(interactions.filepath) == 'character'){
    interactions <- read.interaction.file(interactions.filepath)
  } else if(typeof(interactions.filepath) == 'list'){
    interactions <- read.interaction.files(utr5.filepath = interactions.filepath[[1]],
                                           cds.filepath = interactions.filepath[[2]],
                                           utr3.filepath = interactions.filepath[[3]])
  } else {
    cat(paste("Parameter interactions.filepath is of wrong object type: ", typeof(interactions.filepath), sep = ""))
  }
  sample_annotation <- read.sample.annotation.file(sampleannot.filepath)
  return(list(mir_expr, tpm_expr, interactions, sample_annotation))
}



#-------------------------------------------------------------------------------------
#                               Preprocess methods
#-------------------------------------------------------------------------------------


# #' Search RefSeq mRNA IDs of your interactions data in Ensembl database to obtain matching Ensembl transcript IDs.
# #'
# #' @param interactions Data Frame of miRNA - mRNA interactions. The column name of the mRNAs needs to be "mRNA".
# #' @param ensembl_mart Mart object compatible with the biomaRt query. Gives the possibility to search older versions of the Ensembl database. The default is the current Ensembl database version.
# #' @param id.type Search Ensembl database for "validated", "predicted" or "both" RefSeq mRNA IDs. The default will search for both ID types.
# #' @param log Set a boolean value to print log information to the console. The default is FALSE and the method won't print a log.
# #' @return interactions Modified interactions data with matching Ensembl transcript IDs
preproc.idconversion.refseqtoensembl <- function(interactions, ensembl_mart = NULL, id.type = 'both', log = FALSE){
  
  # Collect input log data
  input.log <- paste("Input interaction contains:\n",
                    "\t- ", length(unique(interactions$mRNA)), " RefSeq mRNA IDs\n", 
                    "\t- ", length(unique(interactions$Genesymbol)), " Genesymbols\n",
                    "\t- ", length(unique(interactions$miRNA)), " miRNA IDs\n",
                    "\t- ", nrow(interactions), " entries\n\n",
                    sep = "")
  
  # Step 1: Prepare reqired data
    refseq_ids <- unique(interactions$mRNA)
    if(is.null(ensembl_mart))
      ensembl_mart <- ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Step 2: Retrieve Ensemble IDs via biomaRt
    mappings <- NULL
    if(id.type == 'validated'){
      mappings <- getBM(attributes = c('refseq_mrna', 'ensembl_transcript_id'),
                        filters = 'refseq_mrna',
                        values = refseq_ids,
                        mart = ensembl_mart)
      colnames(mappings) <- c('mRNA', 'ensembl_transcript_id')
    } else if(id.type == 'predicted'){
      mappings <- getBM(attributes = c('refseq_mrna_predicted', 'ensembl_transcript_id'),
                        filters = 'refseq_mrna_predicted',
                        values = refseq_ids,
                        mart = ensembl_mart)
      colnames(mappings) <- c('mRNA', 'ensembl_transcript_id')
    } else if(id.type == 'both'){
      refseq_ensembl <- getBM(attributes = c('refseq_mrna', 'ensembl_transcript_id'),
                              filters = 'refseq_mrna',
                              values = refseq_ids,
                              mart = ensembl_mart)
      predrefseq_ensembl <- getBM(attributes = c('refseq_mrna_predicted', 'ensembl_transcript_id'),
                                  filters = 'refseq_mrna_predicted',
                                  values = refseq_ids,
                                  mart = ensembl_mart)
      colnames(refseq_ensembl) <- c('mRNA', 'ensembl_transcript_id')
      colnames(predrefseq_ensembl) <- c('mRNA', 'ensembl_transcript_id')
      mappings <- rbind(refseq_ensembl, predrefseq_ensembl)
    } else {
      cat("\nIncorrect value for parameter id_type\nPlease choose between 'validated', 'predicted' and 'both'\nRefSeq to Ensemble IDs conversion stopped...\n")
      return(NULL)
    }
    if(is.null(mappings)){
      cat("\nNo Ensemble IDs found\nRefSeq to Ensemble IDs conversion stopped...\n")
      return(NULL)
    }
  
  # Step 3: Remove entries where no Ensemble IDs were found 
  # and entries where found Genesymbols don't match the interaction data
    blank_interactions <- interactions[!interactions$mRNA %in% mappings$mRNA, ]
    interactions <- interactions[interactions$mRNA %in% mappings$mRNA, ]
    interactions <- merge(interactions, mappings, by = c("mRNA"), allow.cartesian = TRUE)
    
  # LOG
    if(log){
      ambiguous_mappings <- mappings[duplicated(mappings$mRNA), ]
      distinct_interactions <- interactions[!(interactions$mRNA %in% ambiguous_mappings$mRNA), ]
      ambiguous_interactions <- interactions[interactions$mRNA %in% ambiguous_mappings$mRNA, ]
      
      cat("\n*Log of RefSeq ID to Ensembl ID conversion*\n\n")
      
      cat(paste("Used ensembl_mart host: ", listEnsembl(ensembl_mart)$version[1], "\n", sep = ""))
      cat(paste("Used id types: ", id.type, "\n\n", sep = ""))
      
      cat(input.log)
      
      cat("Mapped interactions contain:\n")
      cat(paste("\t- ", length(unique(interactions$mRNA)), " RefSeq mRNA IDs\n", sep = ""))
      cat(paste("\t- ", length(unique(interactions$ensembl_transcript_id)), " Ensembl transcript IDs\n", sep = ""))
      cat(paste("\t- ", length(unique(interactions$Genesymbol)), " Genesymbols\n", sep = ""))
      cat(paste("\t- ", length(unique(interactions$miRNA)), " miRNA IDs\n", sep = ""))
      cat(paste("\t- ", nrow(interactions), " entries\n\n", sep = ""))
      
      cat("Detailed ID mapping log:\n")
      cat(paste("\t- ", length(unique(distinct_interactions$mRNA)), " distictly mapped RefSeq mRNA IDs\n", sep = ""))
      cat(paste("\t- ", length(unique(distinct_interactions$ensembl_transcript_id)), " Ensembl transcript IDs in distinct RefSeq mRNA ID mapping\n", sep = ""))
      cat(paste("\t- ", nrow(distinct_interactions), " entries\n", sep = ""))
      cat(paste("\t- ", length(unique(ambiguous_interactions$mRNA)), " ambiguously mapped RefSeq mRNA IDs\n", sep = ""))
      cat(paste("\t- ", length(unique(ambiguous_interactions$ensembl_transcript_id)), " Ensembl transcript IDs in ambiguous RefSeq mRNA ID mapping\n", sep = ""))
      cat(paste("\t- ", length(unique(ambiguous_interactions$Genesymbol)), " unique Genesymbols in ambiguous RefSeq mRNA ID mapping\n", sep = ""))
      cat(paste("\t- ", nrow(ambiguous_interactions), " entries\n", sep = ""))
      
      cat("ID conversion removed from interactions:\n")
      cat(paste("\t- ", length(unique(blank_interactions$mRNA)), " unique RefSeq mRNA IDs\n", sep = ""))
      cat(paste("\t- ", nrow(blank_interactions), " entries\n", sep = ""))
    }

  return(interactions)
    rm(refseq_ids, mappings, refseq_ensembl, predrefseq_ensembl, blank_interactions, distinct_interactions, ambiguous_interactions, ambiguous_mappings, id.type, log, input.log)
}


preproc.filter.samples <- function(mir_expr, tpm_expr, sample_annotation, primary.disease, sample.type.id){
  # Step 1: Retrieve sample ID list matching primary.disease & sample.type.id
    if(!is.null(sample.type.id)){
      sample.type.id <- as.numeric(sample.type.id)
    }
    
    # Get lists of valid disease types and sample type IDs
    disease_types <- sort(unique(sample_annotation$X_primary_disease))
    sample_types <- unique(subset(sample_annotation, select = c("sample_type_id", "sample_type")))
    sample_types <- sample_types[order(sample_types$sample_type_id), ]
    row.names(sample_types) <- NULL
    
    # Validate input values for primary_disease and sample_type_id
    valid.disease.type <- TRUE
    valid.sample.type <- TRUE
    if(is.null(primary.disease) || missing(primary.disease) || !primary.disease %in% disease_types)
      valid.disease.type <- FALSE
    if(is.null(sample.type.id) || missing(sample.type.id) || !sample.type.id %in% sample_types$sample_type_id)
      valid.sample.type <- FALSE
    
    # Filter tcga sample annotations for primary_disease and sample_type_id
    sample_list <- NULL
    if(!valid.disease.type && !valid.sample.type) {
      cat("\nLIST OF TCGA PRIMARY DISEASE IDENTIFIERS\n")
      cat("----------------------------------------\n\n")
      print(disease_types)
      cat("\n")
      cat("\nLIST OF SAMPLE TYPES WITH INDICES\n")
      cat("---------------------------------\n\n")
      print(sample_types)
      cat("\n")
      return(list(NULL, NULL))
    } else if(!valid.disease.type) {
      sample_list <- sample_annotation[sample_annotation$sample_type_id == sample.type.id & 
                                         !is.na(sample_annotation$sample_type_id), ]
    } else if(!valid.sample.type) {
      sample_list <- sample_annotation[sample_annotation$X_primary_disease == primary.disease, ]
    } else {
      sample_list <- sample_annotation[sample_annotation$X_primary_disease == primary.disease & 
                                         sample_annotation$sample_type_id == sample.type.id & 
                                         !is.na(sample_annotation$sample_type_id), ]
    }
    sample_list <- as.vector(sample_list$sample)
    
    # Report if no samples were found
    if(length(sample_list) < 1){
      cat("\nRetrieved sample list of size 0.\nFilter expr data for samples stopped...\n")
      return(list(NULL, NULL))
    }
  
  # Step 2: Subset mir_expr & gene_expr columns for found sample IDs    
    mir_expr <- mir_expr[ , colnames(mir_expr) %in% sample_list]
    tpm_expr <- tpm_expr[ , colnames(tpm_expr) %in% sample_list]
    
  return(list(mir_expr, tpm_expr))
}

# na.equivalent:  miRNA        0
#                 transcripts  -9.9658
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

preproc.rowsubset.variance <- function(expr_data, top.number = 10){
  expr_data <- as.data.frame(expr_data)
  expr_data$variance <- apply(data.frame(expr_data), 1, var)
  expr_data <- expr_data[order(expr_data$variance, decreasing = T), ]
  expr_data <- expr_data[1:top.number, colnames(expr_data) != 'variance']
  return(expr_data)
}

preproc.colsubset.variance <- function(expr_data, top.number = 10){
  expr_data <- preproc.rowsubset.variance(expr_data = t(tpm_expr),
                                          top.number = top.number)
  expr_data <- t(expr_data)
  return(expr_data)
}

preproc.crossfilter <- function(mir_expr, tpm_expr, interactions){
  mir_expr <- mir_expr[rownames(mir_expr) %in% interactions$miRNA , colnames(mir_expr) %in% colnames(tpm_expr)]
  tpm_expr <- tpm_expr[rownames(tpm_expr) %in% interactions$ensembl_transcript_id, colnames(tpm_expr) %in% colnames(mir_expr)]
  interactions <- interactions[interactions$miRNA %in% rownames(mir_expr) & interactions$ensembl_transcript_id %in% rownames(tpm_expr), ]
  return(list(mir_expr, tpm_expr, interactions))
}

preproc.data.spongefilter <- function(mir_expr, tpm_expr, interactions){
  mir_expr <- t(mir_expr)
  tpm_expr <- t(tpm_expr)
  #rownames(expr_data) <- gsub("\\.", "-", rownames(expr_data))
  interactions_matrix <- as.matrix(as.data.frame.matrix(table(interactions[ , c("ensembl_transcript_id", "miRNA")])))
  return(list(mir_expr, tpm_expr, interactions_matrix))
}


# MASTER METHOD
preproc <- function(mir_expr, tpm_expr, interactions, sample_annotation, primary.disease, sample.type.id, ensembl_mart = NULL, top.number = NULL){
    
  # Step 1: ID conversion  
  interactions <- preproc.idconversion.refseqtoensembl(interactions = interactions,
                                                       ensembl_mart = ensembl_mart)
  
  # Step 2: Extract samples for specific disease and sample type
  box <- preproc.filter.samples(mir_expr = mir_expr,
                                tpm_expr = tpm_expr,
                                sample_annotation = sample_annotation,
                                primary.disease = primary.disease,
                                sample.type.id = sample.type.id)
  
  # Step 3: Filter expression data for NAs and variance
  mir_expr <- preproc.filter.variance(preproc.filter.na(expr_data = box[[1]],
                                                        na.equivalent = 0))
  tpm_expr <- preproc.filter.variance(preproc.filter.na(expr_data = box[[2]], 
                                                        na.equivalent = -9.9658))
  
  # Step 3: Crossfilter data
  box <- preproc.crossfilter(mir_expr = mir_expr,
                             tpm_expr = tpm_expr,
                             interactions = interactions)
  mir_expr <- box[[1]]
  tpm_expr <- box[[2]]
  interactions <- box[[3]]
  
  # Step 4: Subset data if reqired
  if(!is.null(top.number)){
    tpm_expr <- preproc.rowsubset.variance(expr_data = tpm_expr,
                                           top.number = top.number)
    box <- preproc.crossfilter(mir_expr = mir_expr,
                               tpm_expr = tpm_expr,
                               interactions = interactions)
    mir_expr <- box[[1]]
    tpm_expr <- box[[2]]
    interactions <- box[[3]]
      cat(paste("Produced subset of ", top.number, " most variant transcripts\n", sep = ""))
  }
  
  # Step 5: METHOD TANGLING (not usefull but needed for runtime analysis; Redo it!)
  box <- preproc.data.spongefilter(mir_expr = mir_expr,
                                   tpm_expr = tpm_expr,
                                   interactions = interactions)
  mir_expr <- box[[1]]
  tpm_expr <- box[[2]]
  interactions_matrix <- box[[3]]

  return(list(mir_expr, tpm_expr, interactions, interactions_matrix))
}


#-------------------------------------------------------------------------------------
#                               SPONGE methods
#-------------------------------------------------------------------------------------


sponge.unlist.candidates <- function(sponge_filtered){
  candidates <- do.call(rbind, sponge_filtered)
  candidates$transcript <- lapply(strsplit(rownames(candidates), "\\."), '[', 1)
  rownames(candidates) <- NULL
  return(candidates)
}


# MASTER METHOD
sponge.filter <- function(tpm_expr, mir_expr, interactions_matrix, cluster.size = 40){
      start.time <- Sys.time()
      cat("______________________________________________________________________\n")
      cat(paste(start.time, "\t\tStart running the SPONGE interaction filter on ", cluster.size, " cluster...\n", sep = ""))
        cl <- makePSOCKcluster(cluster.size)
        registerDoParallel(cl)

  sponge_filtered_le <- sponge_gene_miRNA_interaction_filter(gene_expr = tpm_expr,
                                                             mir_expr = mir_expr,
                                                             mir_predicted_targets = interactions_matrix,
                                                             coefficient.threshold = -0.05,
                                                             coefficient.direction = "<",
                                                             parallel.chunks = cluster.size)
      cat(paste(Sys.time(), "\t\t\tFinished the SPONGE interaction filter step for negativ correlation\n"), sep = "")
  sponge_filtered_ge <- sponge_gene_miRNA_interaction_filter(gene_expr = tpm_expr,
                                                             mir_expr = mir_expr,
                                                             mir_predicted_targets = interactions_matrix,
                                                             coefficient.threshold = 0.05,
                                                             coefficient.direction = ">",
                                                             parallel.chunks = cluster.size)
      cat(paste(Sys.time(), "\t\t\tFinished the SPONGE interaction filter step for positiv correlation\n"), sep = "\t")
  le_candidates <- sponge.unlist.candidates(sponge_filtered_le)
  le_candidates$correlation <- rep("<", nrow(le_candidates))
  ge_candidates <- sponge.unlist.candidates(sponge_filtered_ge)
  ge_candidates$correlation <- rep(">", nrow(ge_candidates))
  candidates <- rbind(le_candidates, ge_candidates)
  colnames(candidates) <- c('miRNA', 'coefficient', 'ensembl_transcript_id', 'correlation')
  
        stopCluster(cl)
      end.time <- Sys.time()
      diff.time <- difftime(end.time, start.time, units = "secs")
      cat(paste(end.time, "\t\tFinished total run of the SPONGE interaction filter\n=> Total procession time:\t", diff.time, " SECS\n\n", sep = ""))
  
  return(candidates)
}


#-------------------------------------------------------------------------------------
#                                   Position analysis
#-------------------------------------------------------------------------------------

# Finished
exon.model.preproc <- function(candidates, interactions){
  candidates <- merge(candidates, interactions, allow.cartesian = TRUE)
  candidates$binding_site_start <- as.numeric(sapply(strsplit(candidates$binding_site, ","), "[", 1))
  candidates$binding_site_end <- as.numeric(sapply(strsplit(candidates$binding_site, ","), "[", 2))
  candidates$ensembl_transcript_id <- as.character(candidates$ensembl_transcript_id)
  candidates <- select(candidates, miRNA, mRNA, ensembl_transcript_id, Genesymbol, genomic_region, binding_site_start, binding_site_end, binding_probability, coefficient, correlation)
  return(candidates)
}

# Finished
exon.model.build <- function(interactions, ensembl_mart = NULL){
  
  if(is.null(ensembl_mart))
    ensembl_mart <- ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  transcriptomic_annot <- getBM(attributes = c('refseq_mrna', 'strand', 'ensembl_gene_id', 'start_position', 'end_position', 'ensembl_transcript_id', 'ensembl_exon_id', 'transcript_start', 'transcript_end'),
                                filters = 'refseq_mrna',
                                values = unique(interactions$mRNA),
                                mart = ensembl_mart)
  transcriptomic_annot_predicted <- getBM(attributes = c('refseq_mrna_predicted', 'strand', 'ensembl_gene_id', 'start_position', 'end_position', 'ensembl_transcript_id', 'ensembl_exon_id', 'transcript_start', 'transcript_end'),
                                          filters = 'refseq_mrna_predicted',
                                          values = unique(candidates$mRNA),
                                          mart = ensembl_mart)
  colnames(transcriptomic_annot)[1] <- "mRNA"
  colnames(transcriptomic_annot_predicted)[1] <- "mRNA"
  transcriptomic_annot <- rbind(transcriptomic_annot, transcriptomic_annot_predicted)
  
  exonic_annot <- getBM(attributes = c('ensembl_exon_id', 'exon_chrom_start', 'exon_chrom_end'),
                        filters = 'ensembl_exon_id',
                        values = unique(transcriptomic_annot$ensembl_exon_id),
                        mart = ensembl_mart)
  
  annotation <- merge(transcriptomic_annot, exonic_annot, by = "ensembl_exon_id", allow.cartesian = TRUE)
  posstrand <- annotation[annotation$strand == "1", ]
  posstrand$exon_transcript_start <- posstrand$exon_chrom_start - posstrand$transcript_start + 1
  posstrand$exon_transcript_end <- posstrand$exon_chrom_end - posstrand$transcript_start + 1
  posstrand$exon_gene_start <- posstrand$exon_chrom_start - posstrand$start_position + 1
  posstrand$exon_gene_end <- posstrand$exon_chrom_end - posstrand$start_position + 1
  negstrand <- annotation[annotation$strand == "-1", ]
  negstrand$exon_transcript_start <- negstrand$transcript_end - negstrand$exon_chrom_end + 1
  negstrand$exon_transcript_end <- negstrand$transcript_end - negstrand$exon_chrom_start + 1
  negstrand$exon_gene_start <- negstrand$end_position - negstrand$exon_chrom_end + 1
  negstrand$exon_gene_end <- negstrand$end_position - negstrand$exon_chrom_start + 1
  annotation <- rbind(posstrand, negstrand)
  annotation$exon_length <- annotation$exon_transcript_end - annotation$exon_transcript_start + 1
  
  annotation_temp <- NULL
  count <- 0
  pb <- txtProgressBar(min = 0, max = length(unique(annotation$ensembl_transcript_id)), style = 3)
  for(t in unique(annotation$ensembl_transcript_id)){
    sub <- annotation[annotation$ensembl_transcript_id == t, ]
    sub <- sub[order(sub$exon_transcript_start), ]
    sub$exon_mrna_start <- rep(0, nrow(sub))
    sub$exon_mrna_end <- rep(0, nrow(sub))
    index <- 1
    for(i in 1:nrow(sub)){
      sub$exon_mrna_start[i] <- index
      sub$exon_mrna_end[i] <- sub$exon_mrna_start[i] + sub$exon_length[i] - 1
      index <- sub$exon_mrna_end[i] + 1
    }
    if(is.null(annotation_temp))
      annotation_temp <- sub
    else
      annotation_temp <- rbind(annotation_temp, sub)
    
    count <- count + 1
    setTxtProgressBar(pb, count)
  }
  close(pb)
  annotation <- annotation_temp
  annotation <- select(annotation, mRNA, ensembl_gene_id, ensembl_transcript_id, ensembl_exon_id, exon_mrna_start, exon_mrna_end)
  
  return(annotation)
}


# MASTER METHOD
exon.model.enrich <- function(interactions, candidates, ensembl_mart){
  
  # Step 1: Preproc candidates data
  candidates <- exon.model.preproc(candidates = candidates,
                                   interactions = interactions)
  
  # Step 2: Retrieve ensembl annotation data on transcriptomic and exonic positions
  annotation <- exon.model.build(interactions = interactions,
                                  ensembl_mart = ensembl_mart)
  
  # Step 3: Match binding sites to exons
  candidate_exons <- merge(candidates, annotation, by = c("mRNA", "ensembl_transcript_id"), allow.cartesian = TRUE)
  candidate_exons_matched <- candidate_exons %>% filter(exon_mrna_start <= binding_site_start &
                                                        binding_site_end <= exon_mrna_end)
  # => WHAT HAPPENS WITH UNMATCHED BINDING SITES??
  #candidate_exons_unmatched <- anti_join(candidates, candidate_exons_matched, by = c("miRNA", "mRNA", "ensembl_transcript_id", "binding_site_start", "binding_site_end", "correlation"))
  #candidate_exons_unmatched <- merge(candidate_exons_unmatched, annotation, by = c("mRNA", "ensembl_transcript_id"), allow.cartesian = TRUE)
  #candidate_exons_mappable_vers92 <- candidate_exons_unmatched_vers92 %>% filter(exon_mrna_start <= binding_site_start &
  #                                                             binding_site_start <= exon_mrna_end |
  #                                                             exon_mrna_start <= binding_site_end &
  #                                                             binding_site_end <= exon_mrna_end)
  
  # Step 4: Find transcripts without the exons
  mirreti_annotation <- NULL
  for(g in unique(candidate_exons_matched$ensembl_gene_id)){
    sub_annot <- annotation[annotation$ensembl_gene_id == g, ]
    if(length(unique(sub_annot$ensembl_transcript_id)) <= 1)
      next
    sub_candidates <- candidate_exons_matched[candidate_exons_matched$ensembl_gene_id == g, ]
    for(i in 1:nrow(sub_candidates)){
      found_exons <- sub_annot[sub_annot$ensembl_exon_id == sub_candidates$ensembl_exon_id[i], ]
      transcript_annot_with_missing_exon <- sub_annot[!sub_annot$ensembl_transcript_id %in% found_exons$ensembl_transcript_id, ]
      if(length(unique(transcript_annot_with_missing_exon$ensembl_transcript_id)) < 1)
        next
      transcript_annot_with_missing_exon <- transcript_annot_with_missing_exon[transcript_annot_with_missing_exon$ensembl_transcript_id != sub_candidates$ensembl_transcript_id[1], ]
      if(length(unique(transcript_annot_with_missing_exon$ensembl_transcript_id)) > 1){
        dataframe <- data.frame(ensembl_transcript_id_1 = rep(sub_candidates$ensembl_transcript_id[i], length(unique(transcript_annot_with_missing_exon$ensembl_transcript_id))),
                                ensembl_exon_id_1 = rep(sub_candidates$ensembl_exon_id[i], length(unique(transcript_annot_with_missing_exon$ensembl_transcript_id))),
                                correlation = rep(sub_candidates$correlation[i], length(unique(transcript_annot_with_missing_exon$ensembl_transcript_id))))
        dataframe$transcript_without_exon_and_correlation <- unique(transcript_annot_with_missing_exon$ensembl_transcript_id)
        if(is.null(mirreti_annotation)){
          mirreti_annotation <- dataframe
        } else {
          mirreti_annotation <- rbind(mirreti_annotation, dataframe)
        }
      }
    }
  }
  #candidate_exons_matched <- candidate_exons_matched[!candidate_exons_matched$ensembl_gene_id %in% single_transcript_genes, ]
  #annotation <- annotation[!annotation$ensembl_gene_id %in% single_transcript_genes, ]
  return(mirreti_annotation)
}


