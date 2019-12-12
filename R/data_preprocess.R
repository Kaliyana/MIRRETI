library(data.table)
library(graphics)
library(dplyr)
library(biomaRt)


#                                   Parameter
#-------------------------------------------------------------------------------------

primary.disease <- "breast invasive carcinoma"
sample.type.id <- "11"
out.filepath <- "/nfs/home/students/evelyn/bachelor/data/thirdShot/spongerun_data_breastinvasivecarcinome_solidtissuenormal.RDS"



#-------------------------------------------------------------------------------------
#                                 Basic read methods
#-------------------------------------------------------------------------------------

read.expression.file <- function(expression_filepath){
  expr_data <- data.frame(fread(expression_filepath, header = T, sep = "\t"), row.names = 1, check.names = F)
  rownames(expr_data) <- lapply(strsplit(rownames(expr_data), "\\."), '[', 1)
  return(expr_data)
}

read.interaction.file <- function(interaction_filepath){
  interactions <- fread(interaction_filepath, header = T, sep = "\t")
  return(interactions)
}

read.sample.annotation.file <- function(sampleannotation_filepath){
  sample_annotation <- read.csv(sampleannotation_filepath, header = T, sep = "\t")
  return(sample_annotation)
}


#-------------------------------------------------------------------------------------
#                               Preprocess methods
#-------------------------------------------------------------------------------------

# Finished but needs testing
preproc.idconversion.refseqtoensembl(interactions, id.type = 'both'){
  
  # Step 1: Get RefSeq IDs from interactions data
    refseq_ids <- unique(interactions$mRNA)
  
  # Step 2: Retrieve Ensemble IDs via biomaRt
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mappings <- NULL
    if(id.type == 'validated'){
      mappings <- getBM(attributes = c('refseq_mrna', 'ensembl_transcript_id', 'external_gene_name'),
                        filters = 'refseq_mrna',
                        values = refseq_ids,
                        mart = ensembl)
    } else if(id.type == 'predicted'){
      mappings <- getBM(attributes = c('refseq_mrna_predicted', 'ensembl_transcript_id', 'external_gene_name'),
                        filters = 'refseq_mrna_predicted',
                        values = refseq_ids,
                        mart = ensembl)
      colnames(mappings)[1] <- "refseq_mrna"
    } else if(id.type == 'both'){
      refseq_ensembl <- getBM(attributes = c('refseq_mrna', 'ensembl_transcript_id', 'external_gene_name'),
                              filters = 'refseq_mrna',
                              values = refseq_ids,
                              mart = ensembl)
      predrefseq_ensembl <- getBM(attributes = c('refseq_mrna_predicted', 'ensembl_transcript_id', 'external_gene_name'),
                                  filters = 'refseq_mrna_predicted',
                                  values = refseq_ids,
                                  mart = ensembl)
      colnames(refseq_ensembl) <- c('mRNA', 'ensembl_transcript_id', 'Genesymbol')
      colnames(predrefseq_ensembl) <- c('mRNA', 'ensembl_transcript_id', 'Genesymbol')
      mappings <- rbind(refseq_ensembl, predrefseq_ensembl)
    } else {
      cat("\nIncorrect value for parameter id_type\nPlease choose between 'validated', 'predicted' and 'both'\n")
      quit
    }
  
  # Step 3: Remove entries where no Ensemble IDs were found
    mappings <- mappings[!(mappings$ensembl_transcript_id == ""),]
    
  # Step 4: Merge ID mapping to interactions data
    interactions <- interactions[interactions$mRNA %in% mappings$mRNA, ]
    interactions <- merge(interactions, mappings, by = c("mRNA", "Genesymbol"), allow.cartesian = T)
}

# Finished
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
      .Internal(.invokeRestart(list(NULL, NULL), NULL))
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
      cat("\nRetrieved sample list of size 0.\n")
    }
  
  # Step 2: Subset mir_expr & gene_expr columns for found sample IDs    
    mir_expr <- mir_expr[ , colnames(mir_expr) %in% sample_list]
    tpm_expr <- tpm_expr[ , colnames(tpm_expr) %in% sample_list]
}

#Finished
# na.equivalent:  miRNA        0
#                 transcripts  -9.9658
preproc.filter.na <- function(expr_data, na.equivalent, max.na.percent = 0.2){
  na_content <- rowSums(expr_data == na.equivalent)
  expr_data <- expr_data[na_content/ncol(expr_data) <= max.na.percent, ]
  na_content <- colSums(expr_data == na.equivalent)
  expr_data <- expr_data[ , na_content/nrow(expr_data) <= max.na.percent]
}

# Finished
preproc.filter.variance <- function(expr_data, var.threshold = 0.2){
  expr_data$variance <- apply(data.frame(expr_data), 1, var)
  expr_data <- expr_data[order(expr_data$variance, decreasing = T), ]
  threshold <- floor((1-var.threshold) * nrow(expr_data))
  expr_data <- expr_data[1:threshold, colnames(expr_data) != 'variance']
}

# Finished
preproc.subset.variance <- function(expr_data, top.number = 10){
  expr_data$variance <- apply(data.frame(expr_data), 1, var)
  expr_data <- expr_data[order(expr_data$variance, decreasing = T), ]
  expr_data <- expr_data[1:top.number, colnames(expr_data) != 'variance']
}


preproc.crossfilter <- function(mir_expr, tpm_expr, interactions){
  # Step 1: filter for intersecting samples
  
  # Step 2: filter for intersecting miRNAs
  
  # Step 3: filter for intersectiont transcripts
}




#                                    Workflow
#-------------------------------------------------------------------------------------

process.data.for.SPONGEinteractionFilter <- function(mir_expr, tpm_expr, interactions, sample_annotation, primary_disease, sample_type_id, out_filepath){
  #Step 1
  sample_list <- retrieve.tcgasample.list(sample_annotation, primary_disease = primary_disease, sample_type_id = sample_type_id)
  mir_expr <- preprocess.exprdata(mir_expr, 
                                  id_list = unique(interactions$miRNA), 
                                  sample_list = sample_list)
  interactions <- interactions[interactions$miRNA %in% colnames(mir_expr),]
  
  
  #Step 2
  id_mappings <- convert.refseqids.to.ensemblids(unique(interactions$mRNA))
  id_mappings <- id_mappings[id_mappings$ensembl_transcript_id %in% rownames(tpm_expr),]
  
  
  #Step 3
  interactions <- merge(interactions, id_mappings, by = c("mRNA", "Genesymbol"), allow.cartesian = T)
  interactions_matrix <- as.data.frame.matrix(table(select(interactions, "ensembl_transcript_id", "miRNA")))
  tpm_expr <- preprocess.exprdata(tpm_expr,
                                  id_list = unique(interactions$ensembl_transcript_id),
                                  sample_list = sample_list)
  save(as.matrix(mir_expr), as.matrix(tpm_expr), as.matrix(interactions_matrix), file = out_filepath)
}



#-------------------------------------------------------------------------------------
#                                     Methods
#-------------------------------------------------------------------------------------


#' Retrieve a list of TCGA samples filtered for primary disease type and sample type.
#' NOTE: If you don not know which values to use for primary_disease and sample_type_id, use retrieve_TCGA_sample_list(sample_annotation)
#'
#' @param sample_annotation TCGA annotation of samples as data frame.
#' @param primary_disease Name of the cancertype you are interested in.
#' @param sample_type_id A number generated by TCGA as ID for sample types.
#' @return A list of TCGA sample IDs for chosen primary disease type and sample type.
retrieve.tcgasample.list <- function(sample_annotation, primary_disease = NULL, sample_type_id = NULL){

}


preprocess.exprdata <- function(expr_data, id_list, sample_list){
  # Filter expression dataset for molecule IDs in id_list
  rownames(expr_data) <- lapply(strsplit(rownames(expr_data), "\\."), '[', 1)
  expr_data <- expr_data[rownames(expr_data) %in% id_list, ]
  
  # Transpose data set
  expr_data <- as.data.frame(t(expr_data))
  rownames(expr_data) <- gsub("\\.", "-", rownames(expr_data))
  
  # Filter expression dataset for sample IDs in sample_list
  expr_data_filtered <- subset(expr_data, row.names(expr_data) %in% sample_list)
  
  # Spot sample IDs missing in the current expression dataset and 
  # build matrix with values of 0 for the missing sample IDs
  missing_sampleids <- sample_list[!(sample_list %in% row.names(expr_data))]
  missing_expr_data <- as.data.frame(matrix(0,
                                            length(missing_sampleids),
                                            length(colnames(expr_data_filtered))))
  row.names(missing_expr_data) <- missing_sampleids
  colnames(missing_expr_data) <- colnames(expr_data_filtered)
  
  # Append matrix with missing sample IDs to the expression dataset
  extended_expr_data <- rbind(expr_data_filtered, missing_expr_data)
  return(extended_expr_data)
}
