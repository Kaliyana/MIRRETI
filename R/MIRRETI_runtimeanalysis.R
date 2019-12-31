library(rslurm)
library(SPONGE)


#---------------------------------------------------------------------------
#                               Main methode
#---------------------------------------------------------------------------

sponge.runtime <- function(dragonball){
  source("/nfs/home/students/evelyn/bachelor/R_workspace/MIRRETI/R/MIRRETI.R")
  mir.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena"
  tpm.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/tcga_Kallisto_tpm"
  utr5.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_5UTR.txt"
  cds.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_CDS.txt"
  utr3.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_3UTR.txt"
  sampleannot.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/TCGA_phenotype_denseDataOnlyDownload.tsv"
  primary.disease <- "breast invasive carcinoma"
  sample.type.id = "01"
  out.filepath <- "/nfs/home/students/evelyn/bachelor/data/thirdShot/spongefilter_time_table.RDS"
  
  mir_expr <- read.expression.file(mir.filepath)
  tpm_expr <- read.expression.file(tpm.filepath)
  interactions <- read.interaction.files(utr5.filepath = utr5.filepath,
                                         cds.filepath = cds.filepath,
                                         utr3.filepath = utr3.filepath)
  sample_annotation <- read.sample.annotation.file(sampleannot.filepath)
  
  interactions <- preproc.idconversion.refseqtoensembl(interactions)
  box <- preproc.filter.samples(mir_expr, tpm_expr, sample_annotation, primary.disease, sample.type.id)
  mir_expr <- preproc.filter.variance(preproc.filter.na(box[[1]], 0))
  tpm_expr <- preproc.filter.variance(preproc.filter.na(box[[2]], -9.9658))
  box <- preproc.crossfilter(mir_expr, tpm_expr, interactions)
  mir_expr <- box[[1]]
  tpm_expr <- box[[2]]
  interactions <- box[[3]]
  
  # 10, 50, 100, 500, 1000, 1500, 2500, 5000
  sampling <- c(10, 50, 100, 500, 1000, 1500, 2500, 5000)
  runtime_table <- NULL
  
  for(s in sampling){
    t <- preproc.rowsubset.variance(tpm_expr, s)
    box <- preproc.data.spongefilter(mir_expr, t, interactions)
    m <- box[[1]]
    t <- box[[2]]
    i <- box[[3]]
    
    cat(paste("tpm_size:", s, "\n", sep = "\t"))
    cat(paste("tpm_expr:", dim(t), "\n", sep = "\t"))
    
    start.time <- Sys.time()
    candidates <- sponge.filter(t, m, i)
    end.time <- Sys.time()
    df <- data.frame("tpm_size" = s, "runtime" = end.time - start.time)
    
    if(is.null(runtime_table)){
      runtime_table <- df
    } else {
      runtime_table <- rbind(runtime_table, df)
    }
  }
  
  saveRDS(runtime_table, file = out.filepath)
}


#---------------------------------------------------------------------------
#                             Slurm job request
#---------------------------------------------------------------------------

sjob <- slurm_apply(sponge.runtime, params = data.frame("dragonball" = TRUE), jobname = "runtime analysis sponge",
                    nodes = 1, cpus_per_node = 40)
results <- get_slurm_out(sjob, outtype = "raw")
results