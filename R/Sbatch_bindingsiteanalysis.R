library(rslurm)

bindingsite_analysis <- function(dragonball){
  source("/nfs/home/students/evelyn/bachelor/R_workspace/MIRRETI/R/MIRRETI.R")
  load("/nfs/home/students/evelyn/bachelor/data/thirdShot/datapacks/mirreti-4.0_datapack_breastinvasivecarcinoma_01_EnsemblGenes99.RDS")
  interactions <- mirreti.build.correlationtable(candidates = candidates, 
                                                 interactions = interactions)
  saveRDS(interactions, file = "/nfs/home/students/evelyn/bachelor/R_workspace/MIRRETI/sbatch_outfiles/MIRRETI-4.0_interactionsWithSpongeCorrelations_BRCA_01_Ens99.RDS")
}

sjob <- slurm_apply(bindingsite_analysis, params = data.frame(dragonball = TRUE), jobname = "building_correlation_table",
                    nodes = 1, cpus_per_node = 1)
