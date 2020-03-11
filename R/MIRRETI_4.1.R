

coredata <- function(){
  set.seed(01101990)
  mir.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena"
  tpm.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/tcga_Kallisto_tpm"
  interactions.filepaths <- list("/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_5UTR.txt",
                                "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_CDS.txt",
                                "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_3UTR.txt")
  sampleannot.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/TCGA_phenotype_denseDataOnlyDownload.tsv"
  load("~/bachelor/data/core_data/mirreti_coredata.RDS")
  
  primary.diseases <- "breast invasive carcinoma"
  conditions <- "1"
  ensembl_mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "http://jan2020.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")
  ensembl_mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "http://apr2018.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")
  
  rm(mir.filepath, tpm.filepath, interactions.filepaths, sampleannot.filepath)
  rm(primary.disease, sample.type.id)
  rm(ensembl_mart_92, ensembl_mart_99, ensembl_mart)
}


mirreti <- function(tpm_expr, mir_expr, interactions, samples_annot, primary.disease, ensembl_mart, log = FALSE){
      start.time <- 0
      end.time <- 0

  
  # Step 1: Reading the data
      if(log){
        cat("\n\t\t\t\t*** Welcome to MIRRETI 4.1 ***\n\n")
        start.time <- Sys.time()
        cat("______________________________________________________________________\n")
        cat(paste(start.time, "\t\tStart reading the data...\n", sep = ""))
      }
  data <- mirreti.readData(tpm.filepath = tpm.filepath,
                           mir.filepath = mir.filepath,
                           interactions.filepath = interactions.filepath,
                           sampleannot.filepath = sampleannot.filepath,
                           log = log)
  tpm_expr <- data[[1]]
  mir_expr <- data[[2]]
  interactions <- data[[3]]
  sample_annotation <- data[[4]]
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\t\tSUCCESSFUL\n", sep = ""))
        cat(paste("=> Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
      }

}

mirreti.readData <- function(tpm.filepath, mir.filepath, interactions.filepath, sampleannot.filepath, log = FALSE){
  source(paste(getwd(), "/R/mirreti.readDataformats.R", sep = ""))
      start.time <- 0
      end.time <- 0

  
  # Step 1: read TPM expression data
      if(log){
        cat("\n\nMIRRETI data reading report\n")
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart reading TPM expression data...\n", sep = ""))
      }
  tpm_expr <- read.expressionFile(tpm.filepath)
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\tSUCCESSFUL\n", sep = ""))
        cat(paste("Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
      }
  
  
  # Step 2: read miRNA expression data
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart reading miRNA expression data...\n", sep = ""))
      }
  mir_expr <- read.expressionFile(mir.filepath)
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\tSUCCESSFUL\n", sep = ""))
        cat(paste("Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
      }
  
  
  # Step 3: read mRNA-miRNA interactions data
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart reading mRNA-miRNA interactions data...\n", sep = ""))
      }
  interactions <- NULL
  if(typeof(interactions.filepath) == 'character'){
    interactions <- read.interactionFile(interactions.filepath)
  } else if(typeof(interactions.filepath) == 'list'){
    interactions <- read.interactionFiles(utr5.filepath = interactions.filepath[[1]],
                                          cds.filepath = interactions.filepath[[2]],
                                          utr3.filepath = interactions.filepath[[3]])
  } else {
    cat(paste("Parameter interactions.filepath is of wrong object type: ", typeof(interactions.filepath), sep = ""))
  }
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\tSUCCESSFUL\n", sep = ""))
        cat(paste("Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
      }
  
  
  # Step 4: read TCGA sample annotation
      if(log){
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart reading TCGA sample annotation...\n", sep = ""))
      }
  sample_annotation <- read.sampleAnnotationFile(sampleannot.filepath)
      if(log){
        end.time <- Sys.time()
        cat(paste(end.time, "\tSUCCESSFUL\n", sep = ""))
        cat(paste("Time required:\t", difftime(end.time, start.time, units = "secs"), " secs\n\n", sep = ""))
        cat("End of the MIRRETI data reading report\n\n")
      }
  
  
  return(list(tpm_expr, mir_expr, interactions, sample_annotation))
}

mirreti.preproc <- function(tpm_expr, mir_expr, interactions, sample_annotation, primary.disease, ensembl_mart, log = FALSE){
  source(paste(getwd(), "/R/mirreti.preproc.R", sep = ""))
      start.time <- 0
      end.time <- 0
  
  # Step 1: Filter expression data and TCGA samples for cancer type and condition 01 and 11
      if(log){
        cat("\n\nMIRRETI data preprocessing report\n")
        start.time <- Sys.time()
        cat(paste(start.time, "\tStart filtering expression data and TCGA samples for cancer type: ", primary.disease, " and conditions 01/11...\n", sep = ""))
      }
  data <- preproc.filter.samples(tpm_expr = tpm_expr,
                                 mir_expr = mir_expr,
                                 sample_annotation = sample_annotation,
                                 primary.disease = primary.disease,
                                 conditions = conditions,
                                 log = TRUE)
    
}

mirreti.drimseq <- function(){
  
}

mirreti.spongeInteractionFilter(){
  
}