library(rslurm)
cluster.size = 25

sbatch.mirreti.twoConditions <- function(primary.diseases, conditions, cluster.size){
   set.seed(01101990)
   source("/nfs/home/students/evelyn/bachelor/R_workspace/MIRRETI/R/mirreti.complexrun.R")
   load("~/bachelor/data/core_data/mirreti_coredata.RDS")
   
   data <- tcgaSamples.subset(tpm_expr = tpm_expr,
                              mir_expr = mir_expr,
                              sample_annotation = sample_annotation,
                              primary.diseases = primary.diseases,
                              conditions = conditions,
                              log = TRUE)
   tpm_expr <- data[[1]]
   mir_expr <- data[[2]]
   sample_annotation <- data[[3]]
     rm(data)
     gc()
   interactions <- mirwalkInteractions.idconversion(interactions = interactions)
   
   data <- mirretiData.filter.actualTargets(tpm_expr = tpm_expr,
                                            mir_expr = mir_expr,
                                            interactions = interactions)
   tpm_expr <- data[[1]]
   mir_expr <- data[[2]]
   interactions <- data[[3]]
     rm(data)
     gc()
   
   data <- mirreti.twoConditions(tpm_expr = tpm_expr,
                                 mir_expr = mir_expr,
                                 sample_annotation = sample_annotation,
                                 interactions = interactions,
                                 cluster.size = cluster.size,
                                 log = TRUE)
   interactions <- data[[1]]
   drimseq <- data[[2]]
     rm(data)
     gc()
   
   outfile <- paste("/nfs/home/students/evelyn/bachelor/data/datapacks/mirreti-4.1-_datapack_", paste(gsub(" ", "", primary.diseases), collapse = "|"), "_", paste(conditions, collapse = "|"), ".RDS", sep = "")
   save(tpm_expr, mir_expr, sample_annotation, interactions, drimseq, file = outfile)
   outfile <- paste("/nfs/home/students/evelyn/bachelor/data/datapacks/mirreti-4.1-_minipack_", paste(gsub(" ", "", primary.diseases), collapse = "|"), "_", paste(conditions, collapse = "|"), ".RDS", sep = "")
   save(interactions, drimseq, file = outfile)
}

params <- data.frame(primary.diseases = "breast invasive carcinoma", 
                     conditions = c("1", "11"),
                     cluster.size = cluster.size)

sjob <- slurm_apply(sbatch.mirreti.twoConditions, params = params, jobname = "mirreti_4.1_BRCA_1|11",
                    nodes = 1, cpus_per_node = cluster.size)