library(rslurm)
cluster.size = 25

sbatch.mirreti.simplerun <- function(primary.diseases, conditions, cluster.size){
  set.seed(01101990)
  source("/nfs/home/students/evelyn/bachelor/R_workspace/MIRRETI/R/mirreti.simplerun.R")
  load("~/bachelor/data/core_data/mirreti_coredata.RDS")
  
  data <- tcgaSamples.subset(tpm_expr = tpm_expr,
                             mir_expr = mir_expr,
                             sample_annotation = sample_annotation,
                             primary.diseases = primary.diseases,
                             conditions = conditions,
                             log = TRUE)
  tpm_expr <- data[[1]]
  mir_expr <- data[[2]]
  interactions <- mirwalkInteractions.idconversion(interactions = interactions)
  
  data <- mirretiData.filter.actualTargets(tpm_expr = tpm_expr,
                                           mir_expr = mir_expr,
                                           interactions = interactions)
  tpm_expr <- data[[1]]
  mir_expr <- data[[2]]
  interactions <- data[[3]]
    rm(data)
    gc()
    
  interactions <- mirreti.simplerun(tpm_expr = tpm_expr,
                                    mir_expr = mir_expr,
                                    interactions = interactions,
                                    cluster.size = cluster.size,
                                    log = TRUE)
  
  outfile <- paste("/nfs/home/students/evelyn/bachelor/data/datapacks/mirreti-4.1-simplerun_datapack_", gsub(" ", "", primary.diseases), "_", conditions, ".RDS", sep = "")
  save(tpm_expr, mir_expr, interactions, file = outfile)
}
params <- data.frame(primary.diseases = "breast invasive carcinoma", 
                     conditions = "1",
                     cluster.size = cluster.size)

sjob <- slurm_apply(sbatch.mirreti.simplerun, params = params, jobname = "mirreti_4.1_simplerun",
                    nodes = 1, cpus_per_node = cluster.size)
