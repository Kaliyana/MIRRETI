library(rslurm)
library(SPONGE)
library(doParallel)
library(ggplot2)
library(dplyr)

#sponge.runtime(40)
#spongefilter_runtimes <-  add.to.runtimes()

# HARDCODED!
add.to.runtimes <- function(){
  spongefilter_timetable_temp <- readRDS("~/bachelor/data/thirdShot/spongefilter_timetable_temp.RDS")
  spongefilter_runtimes <- rbind(spongefilter_runtimes, spongefilter_timetable_temp)
  return(spongefilter_runtimes)
}


# HARDCODED!
plot.runtimes <- function(runtimes){
  cluster.size = 40
  runtimes <- spongefilter_runtimes
  runtimes$runtime <- as.numeric(runtimes$runtime)
  runtimes <- runtimes[runtimes$clusters == cluster.size & (runtimes$runtime < 441 | runtimes$runtime > 442), ]
  lm <- lm(runtime ~ tpm_size, data = runtimes)
  nlm <- lm(runtime ~ tpm_size + I(tpm_size*log(tpm_size)), data=runtimes)
  runtimes <- mutate(runtimes, model = predict(nlm))
  ggplot(runtimes, aes(x = tpm_size, y = runtime)) + geom_point() + geom_smooth(method = "lm", se = F) + ggtitle(paste('SPONGE runtime analysis on ', cluster.size, ' cluster - lm fitting', sep = ""))
  ggplot(runtimes) + geom_point(aes(tpm_size, runtime)) + geom_line(aes(tpm_size,model), color="blue")  + ggtitle(paste('SPONGE runtime analysis on ', cluster.size, ' cluster - nlm fitting', sep = ""))
  # tpm_size <- 15834
  tpm_size <- 10000
  predict(lm, newdata = list(tpm_size=tpm_size))
  predict(nlm, newdata = list(tpm_size=tpm_size))
}

#---------------------------------------------------------------------------
#                               Main methode
#---------------------------------------------------------------------------

# OUTDATED!
sponge.runtime <- function(cluster.size){
    cat("\n\t*** Welcome to the MIRRETI runtime analysis ***\n\n")
    source("/nfs/home/students/evelyn/bachelor/R_workspace/MIRRETI/R/MIRRETI.R")
  mir.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena"
  tpm.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/tcga_Kallisto_tpm"
  utr5.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_5UTR.txt"
  cds.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_CDS.txt"
  utr3.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_3UTR.txt"
  sampleannot.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/TCGA_phenotype_denseDataOnlyDownload.tsv"
  primary.disease <- "breast invasive carcinoma"
  sample.type.id = "01"
  out.filepath <- "/nfs/home/students/evelyn/bachelor/data/thirdShot/spongefilter_timetable_temp.RDS"
  cluster.size <- 10
  
    start.time <- Sys.time()
    cat("______________________________________________________________________\n")
    cat(paste(start.time, "Start reading the data...\n", sep = "\t"))
  mir_expr <- read.expression.file(mir.filepath)
  tpm_expr <- read.expression.file(tpm.filepath)
  interactions <- read.interaction.files(utr5.filepath = utr5.filepath,
                                         cds.filepath = cds.filepath,
                                         utr3.filepath = utr3.filepath)
  sample_annotation <- read.sample.annotation.file(sampleannot.filepath)
    end.time <- Sys.time()
    diff.time <- difftime(end.time, start.time, units = "secs")
    cat(paste(end.time, "\tFinished reading the data\nTotal procession time:\t", diff.time, " SECS\n\n", sep = ""))
  
  
    start.time <- Sys.time()
    cat("______________________________________________________________________\n")
    cat(paste(start.time, "Start preprocessing the data...\n", sep = "\t"))
  interactions <- preproc.idconversion.refseqtoensembl(interactions)
  box <- preproc.filter.samples(mir_expr, tpm_expr, sample_annotation, primary.disease, sample.type.id)
  mir_expr <- preproc.filter.variance(preproc.filter.na(box[[1]], 0))
  tpm_expr <- preproc.filter.variance(preproc.filter.na(box[[2]], -9.9658))
  box <- preproc.crossfilter(mir_expr, tpm_expr, interactions)
  mir_expr <- box[[1]]
  tpm_expr <- box[[2]]
  interactions <- box[[3]]
    end.time <- Sys.time()
    diff.time <- difftime(end.time, start.time, units = "secs")
    cat(paste(end.time, "\tFinished preprocessing the data\nTotal procession time:\t", diff.time, " SECS\n\n", sep = ""))
  
  
  # 10, 50, 100, 500, 1000, 1500, 2500, 5000
  sampling <- c(10)
  runtime_table <- NULL
  
  for(s in sampling){
      start.time <- Sys.time()
      cat("______________________________________________________________________\n")
      cat(paste(start.time, "\tStart subsetting the data for ", s, " most variant transcripts...\n", sep = ""))
    t <- preproc.rowsubset.variance(tpm_expr, s)
    box <- preproc.data.spongefilter(mir_expr, t, interactions)
    m <- box[[1]]
    t <- box[[2]]
    i <- box[[3]]
      end.time <- Sys.time()
      diff.time <- difftime(end.time, start.time, units = "secs")
      cat(paste(end.time, "\tFinished subsetting the data for ", s, " most variant transcripts\nTotal procession time:\t", diff.time, " SECS\n", sep = ""))
      cat(paste(" tpm_size:", s, "\n", sep = "\t"))
      cat(paste(" tpm_expr:", dim(t), "\n", sep = "\t"))
      cat("\n")
    
    
      start.time <- Sys.time()
      cat("______________________________________________________________________\n")
      cat(paste(start.time, "\tStart sponge filter for ", s, " most variant transcripts on ", cluster.size, " cluster...\n", sep = ""))
        cl <- makePSOCKcluster(cluster.size)
        registerDoParallel(cl)
    candidates <- sponge.filter(t, m, i, cluster.size)
        stopCluster(cl)
      end.time <- Sys.time()
      diff.time <- difftime(end.time, start.time, units = "secs")
      cat(paste(end.time, "\tFinished sponge filter for ", s, " most variant transcripts\nTotal procession time:\t", diff.time, " SECS\n\n", sep = ""))
    df <- data.frame("tpm_size" = s, "runtime" = diff.time, "clusters" = cluster.size)
    
    if(is.null(runtime_table)){
      runtime_table <- df
    } else {
      runtime_table <- rbind(runtime_table, df)
    }
  }
  
  saveRDS(runtime_table, file = out.filepath)
  cat("\t*** END OF THE MIRRETI RUNTIME ANALYSIS ***")
}


#---------------------------------------------------------------------------
#                             Slurm job request
#---------------------------------------------------------------------------
sjob <- slurm_apply(sponge.runtime, params = data.frame("cluster.size" = 25), jobname = "runtime_analysis_sponge_6",
                    nodes = 1, cpus_per_node = 25)
results <- get_slurm_out(sjob, outtype = "raw")
results