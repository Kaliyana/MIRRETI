library(data.table)   # MANDATORY! - fread()
library(graphics)
library(tidyverse)
library(dplyr)        # MANDATORY! - select(), mutate_all()
library(biomaRt)      # MANDATORY! - useMart(), getBM(), listMarts()
library(SPONGE)       # MANDATORY!
library(doParallel)   # MANDATORY!
library(plyr)


coredata <- function(){
  set.seed(01101990)
  mir.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena"
  tpm.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/tcga_Kallisto_tpm"
  interactions.filepath <- list("/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_5UTR.txt",
                                 "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_CDS.txt",
                                 "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_3UTR.txt")
  sampleannot.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/TCGA_phenotype_denseDataOnlyDownload.tsv"
  primary.disease <- "breast invasive carcinoma"
  sample.type.id <- "01"
  ensembl_mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "http://jan2020.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")
  ensembl_mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "http://apr2018.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")
  
  rm(mir.filepath, tpm.filepath, interactions.filepath, sampleannot.filepath)
  rm(primary.disease, sample.type.id)
  rm(ensembl_mart_92, ensembl_mart_99, ensembl_mart)
}


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
  data <- preproc(mir_expr = show_mir_expr,
                  tpm_expr = show_tpm_expr,
                  interactions = show_interactions,
                  sample_annotation = sample_annotation,
                  primary.disease = primary.disease,
                  sample.type.id = sample.type.id,
                  ensembl_mart = ensembl_mart, 
                  top.number = NULL,
                  log = TRUE)
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
                              interactions_matrix =  interactions_matrix,
                              cluster.size = cluster.size)

      
      cat("\t*** END OF THE MIRRETI RUNTIME ANALYSIS ***\n\n")
  
  return(list(candidates, interactions, mir_expr, tpm_expr))
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
# #' @return interactions Modified interactions data with matching Ensembl transcript IDs
preproc.idconversion.refseqtoensembl <- function(interactions, ensembl_mart = NULL, id.type = 'both'){
  
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
    ambiguous_mrnas <- mappings[duplicated(mappings$mRNA), ]
    ambiguous_transcripts <- mappings[duplicated(mappings$ensembl_transcript_id), ]
    distinct_mappings <- mappings[!(mappings$mRNA %in% ambiguous_mrnas$mRNA) &
                                  !(mappings$ensembl_transcript_id %in% ambiguous_transcripts$ensembl_transcript_id), ]
  # This MIRRETI VERSION USES ONLY DISTINCT MAPPINGS  
    distinct_interactions <- interactions[interactions$mRNA %in% distinct_mappings$mRNA, ]
    interactions <- merge(distinct_interactions, distinct_mappings, by = c("mRNA"), allow.cartesian = TRUE)

  return(interactions)
    rm(refseq_ids, mappings, refseq_ensembl, predrefseq_ensembl, ambiguous_mrnas, ambiguous_transcripts, distinct_interactions, distinct_mappings, id.type)
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
preproc <- function(mir_expr, tpm_expr, interactions, sample_annotation, primary.disease, sample.type.id, ensembl_mart = NULL, top.number = NULL, log = FALSE){
    
  if(log){
    cat("\nMIRRETI data preprocess report\n")
    cat("-Summary of the input data-\n")
    cat(paste("\tmir_expr: ", nrow(mir_expr), " miRNAs and ", ncol(mir_expr), " samples\n"), sep = "")
    cat(paste("\ttpm_expr: ", nrow(tpm_expr), " transcripts and ", ncol(tpm_expr), " samples\n"), sep = "")
    cat(paste("\tinteractions: ", uniqueN(interactions$miRNA), " miRNAs, ", uniqueN(interactions$mRNA), " mRNAs and ", nrow(interactions), "miRNA-mRNA interactions\n"), sep = "")
    cat(paste("\tprimary disease: ", primary.disease, "\n"), sep = "")
    cat(paste("\tsample condition: ", sample.type.id, "\n"), sep = "")
    cat(paste("\tEnsembl version: ", listMarts(ensembl_mart)[1,2], "\n\n"), sep = "")
  }
  
  # Step 1: ID conversion  
  interactions <- preproc.idconversion.refseqtoensembl(interactions = interactions,
                                                       ensembl_mart = ensembl_mart)
  if(log){
    cat("-Log of RefSeq to Ensembl ID conversion of the interactions data-\n")
    cat(paste("\tmiRNAs: ", uniqueN(interactions$miRNA), "\n"), sep = "")
    cat(paste("\tmRNAs: ", uniqueN(interactions$mRNA), "\n"), sep = "")
    cat(paste("\ttranscripts: ", uniqueN(interactions$ensembl_transcript_id), "\n"), sep = "")
    cat(paste("\tmiRNA-mRNA interactions: ", nrow(interactions), "\n\n"), sep = "")
  }
  
  # Step 2: Extract samples for specific disease and sample type
  box <- preproc.filter.samples(mir_expr = mir_expr,
                                tpm_expr = tpm_expr,
                                sample_annotation = sample_annotation,
                                primary.disease = primary.disease,
                                sample.type.id = sample.type.id)
  if(log){
    cat("-Log of expression data sample extraction based on primary disease and sample condition-\n")
    cat(paste("\tmir_expr: ", nrow(box[[1]]), " miRNAs and ", ncol(box[[1]]), " samples\n"), sep = "")
    cat(paste("\ttpm_expr: ", nrow(box[[2]]), " transcripts and ", ncol(box[[2]]), " samples\n\n"), sep = "")
  }
  
  # Step 3: Filter expression data for NAs and variance
  mir_expr <- preproc.filter.variance(preproc.filter.na(expr_data = box[[1]],
                                                        na.equivalent = 0))
  tpm_expr <- preproc.filter.variance(preproc.filter.na(expr_data = box[[2]], 
                                                        na.equivalent = -9.9658))
  if(log){
    cat("-Log of filtering expression data for NAs and variance-\n")
    cat(paste("\tmir_expr: ", nrow(mir_expr), " miRNAs and ", ncol(mir_expr), " samples\n"), sep = "")
    cat(paste("\ttpm_expr: ", nrow(tpm_expr), " transcripts and ", ncol(tpm_expr), " samples\n\n"), sep = "")
  }
  
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
  if(log){
    cat("Summary of the output data-\n")
    cat(paste("\tmir_expr: ", ncol(mir_expr), " miRNAs and ", nrow(mir_expr), " samples\n"), sep = "")
    cat(paste("\ttpm_expr: ", ncol(tpm_expr), " transcripts and ", nrow(tpm_expr), " samples\n"), sep = "")
    cat(paste("\tinteractions: ", uniqueN(interactions$miRNA), " miRNAs, ", uniqueN(interactions$mRNA), " mRNAs, ", uniqueN(interactions$ensembl_transcript_id), " transcripts and ", nrow(interactions), "miRNA-mRNA interactions\n\n"), sep = "")
  }
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
#                                   
#-------------------------------------------------------------------------------------


mirreti.build.correlationtable <- function(candidates, interactions){
  interactions_mrnamirna_pairs <- as.character(unique(paste(interactions$miRNA, interactions$ensembl_transcript_id, sep = ",")))
  interactions_correlated <- as.character(unique(paste(candidates$miRNA, candidates$ensembl_transcript_id, sep = ",")))
  #interactions_correlated <- interactions_correlated %>% mutate_all(as.character)
  interactions_notCorrelated <- interactions_mrnamirna_pairs[!interactions_mrnamirna_pairs %in% interactions_correlated]
  interactions_notCorrelated <- data.frame("miRNA" = lapply(strsplit(interactions_notCorrelated, ','), '[', 1),
                                           "ensembl_transcript_id" = lapply(strsplit(interactions_notCorrelated, ','), '[', 2),
                                           "coefficient" = rep(0, length(interactions_notCorrelated)),
                                           "correlation" = rep("-", length(interactions_notCorrelated)))
  interactions_notCorrelated$miRNA <- as.character(interactions_notCorrelated$miRNA)
  interactions_notCorrelated$ensembl_transcript_id <- as.character(interactions_notCorrelated$ensembl_transcript_id)
  correlation_table <- rbind(candidates, interactions_notCorrelated)
  correlation_table$miRNA <- as.character(correlation_table$miRNA)
  correlation_table$ensembl_transcript_id <- as.character(correlation_table$ensembl_transcript_id)
  correlation_table <- merge(interactions, correlation_table, by = c("miRNA", "ensembl_transcript_id"))
  return(correlation_table)
      rm(interactions_mrnamirna_pairs, interactions_correlated, interactions_notCorrelated)
}

# find genes which lose binding sites within isoforms
mirreti.filter.correlationtable <- function(correlation_table){
  correlation_table_filtered <- data.frame(matrix(ncol = ncol(correlation_table), nrow = 0))
  colnames(correlation_table_filtered) <- colnames(correlation_table)
  
  for(g in unique(correlation_table$Genesymbol)){
    sub <- correlation_table[correlation_table$Genesymbol == g, ]
    
    if(uniqueN(g_sub$ensembl_transcript_id) <= 1){
      next()
    }
    
    correlation_table_filtered <- rbind(correlation_table_filtered, sub)
  }
  
  return(correlation_table_filtered)
}



tryouts <- function(){
  psi_table <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(psi_table) <- c("Genesymbol", "miRNA", "correlation", "psi", "transcripts", "variants")
  
  
  
  
  
  for(g in unique(correlation_table$Genesymbol)){
    g_sub <- correlation_table[correlation_table$Genesymbol == g, ]
    
    unique(g_sub$correlation)
    unique(g_sub$ensembl_transcript_id)
    t1 <- g_sub[g_sub$ensembl_transcript_id == "ENST00000006658", ]
    t2 <- g_sub[g_sub$ensembl_transcript_id == "ENST00000356488", ]
    t5 <- g_sub[g_sub$ensembl_transcript_id == "", ]
    t3 <- t1[!t1$miRNA %in% t2$miRNA, ]
    t4 <- t2[!t2$miRNA %in% t1$miRNA, ]
    
    if(uniqueN(g_sub$ensembl_transcript_id) <= 1 | unique(g_sub$correlation) == "-")
      next()
    for(m in unique(g_sub$miRNA)){
      m_sub <- g_sub[g_sub$miRNA == m, ]
      
      for(c in unique(m_sub$correlation)){
        df <- data.frame("Genesymbol" = m_sub$Genesymbol[1],
                         "miRNA" = m_sub$miRNA[1],
                         "correlation" = unique(m_sub$correlation))
      }
      
    }
    correlation_table_multitranscript <- rbind(correlation_table_multitranscript, sub)
  }
}











#-------------------------------------------------------------------------------------
#                                   Position analysis
#-------------------------------------------------------------------------------------


suppa.preprocess <- function(tpm_expr, gtf, out.filepath){
  interaction_transcripts <- unique(tpm_expr$ensembl_transcript_id)
  gtf <- gtf[gtf$feature == "exon", ]
  temp <- separate(g, attribute, c("gene_id", "gene_version", "transcript_id", "transcript_version", "exon_number", "gene_name", "gene_source", "gene_biotype", "transcript_name", "transcript_source", "transcript_biotype", "exon_id", "exon_version", "tag", "transcript_support_level"), sep = ";")
  temp$transcript_id <- as.character(lapply(strsplit(temp$transcript_id, ' |\\.|"'), '[', 4))
  g <- gtf[temp$transcript_id %in% interaction_transcripts, ]
}

# SUPER COOL METHOD BUT NOT TO BE USED ANYMORE
buildgtf <- function(ensembl_transcript_ids, ensembl_mart = NULL){
  transcriptomic_gtf <- getBM(attributes = c('chromosome_name', 'source', 'strand', 'ensembl_gene_id_version', 'ensembl_transcript_id_version', 'external_gene_name', 'gene_biotype', 'external_transcript_name', 'transcript_source', 'transcript_biotype', 'ccds', 'ensembl_exon_id', 'transcript_tsl'),
               filters = 'ensembl_transcript_id',
               values = unique(ensembl_transcript_ids),
               mart = ensembl_mart)
  exonic_gtf <- getBM(attributes = c('ensembl_exon_id', 'exon_chrom_start', 'exon_chrom_end', 'rank', 'ensembl_transcript_id_version'),
                      filters = 'ensembl_exon_id',
                      values = unique(transcriptomic_gtf$ensembl_exon_id),
                      mart = ensembl_mart)
  
  gtf <- merge(transcriptomic_gtf, exonic_gtf, by = c('ensembl_transcript_id_version', 'ensembl_exon_id'))
  gtf$gene_version <- lapply(strsplit(transcriptomic_gtf$ensembl_gene_id_version, '\\.'), '[', 2)
  gtf$ensembl_gene_id_version <- lapply(strsplit(transcriptomic_gtf$ensembl_gene_id_version, '\\.'), '[', 1)
  gtf$transcript_version <- lapply(strsplit(transcriptomic_gtf$ensembl_transcript_id_version, '\\.'), '[', 2)
  gtf$ensembl_transcript_id_version <- lapply(strsplit(transcriptomic_gtf$ensembl_transcript_id_version, '\\.'), '[', 1)
  gtf$strand <- revalue(as.character(transcriptomic_gtf$strand), c('1' = '+'))
  gtf$strand <- revalue(as.character(transcriptomic_gtf$strand), c('-1' = '-'))
  gtf$feature <- rep('exon', nrow(transcriptomic_gtf))
  gtf$score <- rep('.', nrow(transcriptomic_gtf))
  gtf$frame <- rep('.', nrow(transcriptomic_gtf))
  gtf$transcript_tsl <- gsub("tsl", "", gtf$transcript_tsl)
  gtf$attribute <- paste('gene_id "', gtf$ensembl_gene_id_version, '"; gene_version "', gtf$gene_version, '"; transcript_id "', gtf$ensembl_transcript_id_version, '"; transcript_version "', gtf$transcript_version,
                         '"; gene_name "', gtf$external_gene_name, '"; gene_source "', gtf$source, '"; gene_biotype "', gtf$gene_biotype, '"; transcript_name "', gtf$external_transcript_name, '"; transcript_source "',
                         gtf$transcript_source, '"; transcript_biotype "', gtf$transcript_biotype, '"; ccds_id "', gtf$ccds, '"; exon_id "', gtf$ensembl_exon_id, '"; exon_version "', gtf$rank,
                         '"; transcript_support_level "', gtf$transcript_tsl, '";', sep = "")

  gtf <- select(gtf, chromosome_name, source, feature, exon_chrom_start, exon_chrom_end, score, strand, frame, attribute)
  colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  # out.filepath <- "/nfs/home/students/evelyn/bachelor/data/thirdShot/suppa/GRCh38.99_ensembl_breastcancer_01/GRCh38.99_ensembl_breastcancer_01"
  write.table(gtf, file = paste(out.filepath, ".gtf", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(t(tpm_expr), file = paste(out.filepath, ".tpm", sep = ""), sep = "\t", quote = FALSE)
  }


#-------------------------------------------------------------------------------------
#                                   Position analysis
#-------------------------------------------------------------------------------------

# Finished
bsmapping.preproc <- function(candidates, interactions){
  candidates <- merge(candidates, interactions, allow.cartesian = TRUE)
  candidates$binding_site_start <- as.numeric(sapply(strsplit(candidates$binding_site, ","), "[", 1))
  candidates$binding_site_end <- as.numeric(sapply(strsplit(candidates$binding_site, ","), "[", 2))
  candidates$ensembl_transcript_id <- as.character(candidates$ensembl_transcript_id)
  candidates <- select(candidates, miRNA, mRNA, ensembl_transcript_id, Genesymbol, genomic_region, binding_site_start, binding_site_end, binding_probability, coefficient, correlation)
  return(candidates)
}

#MASTER METHOD
bsmapping.ensembl <- function(candidates, interactions, ensembl_mart = NULL){
  
  candidates <- bsmapping.preproc(candidates = candidates,
                                  interactions = interactions)
  
  if(is.null(ensembl_mart))
    ensembl_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
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
  
  candidate_exons <- merge(candidates, annotation, by = c("mRNA", "ensembl_transcript_id"), allow.cartesian = TRUE)
  candidate_exons_matched <- candidate_exons %>% filter(exon_mrna_start <= binding_site_start &
                                                          binding_site_end <= exon_mrna_end)
  return(candidate_exons_matched)
    rm(transcriptomic_annot, transcriptomic_annot_predicted, exonic_annot)
    rm(posstrand, negstrand, count, pb, i, index, t, sub, annotation_temp)
    rm(annotation, candidate_exons)
}


# DEAD METHOD
exon.model.enrich <- function(interactions, candidates, ensembl_mart){
  
  # Step 1: 
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
        transcripts.without.correlation <- length(unique(transcript_annot_with_missing_exon$ensembl_transcript_id))
        dataframe <- data.frame(mRNA = rep(sub_candidates$mRNA[i], transcripts.without.correlation),
                                miRNA = rep(sub_candidates$miRNA[i], transcripts.without.correlation),
                                ensembl_transcript_id_1 = rep(sub_candidates$ensembl_transcript_id[i], transcripts.without.correlation),
                                ensembl_exon_id_1 = rep(sub_candidates$ensembl_exon_id[i], transcripts.without.correlation),
                                correlation = rep(sub_candidates$correlation[i], transcripts.without.correlation))
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


