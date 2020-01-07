library(rslurm)


#---------------------------------------------------------------------------
#                               Main methode
#---------------------------------------------------------------------------


spongefilter <- function(dragonball){
  source("/nfs/home/students/evelyn/bachelor/R_workspace/MIRRETI/R/MIRRETI.R")
  mir.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena"
  tpm.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/tcga_Kallisto_tpm"
  utr5.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_5UTR.txt"
  cds.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_CDS.txt"
  utr3.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_3UTR.txt"
  sampleannot.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/TCGA_phenotype_denseDataOnlyDownload.tsv"
  primary.disease <- "breast invasive carcinoma"
  sample.type.id = "01"

  candidates <- mirreti(mir.filepath, tpm.filepath, list(utr5.filepath, cds.filepath, utr3.filepath),
                        sampleannot.filepath, primary.disease, sample.type.id, cluster.size = 25, top.number = 100)
  
  saveRDS(candidates, file = "/nfs/home/students/evelyn/bachelor/data/thirdShot/SPONGE_candidates_2.RDS")
}



#---------------------------------------------------------------------------
#                             Slurm job request
#---------------------------------------------------------------------------


sjob <- slurm_apply(spongefilter, params = data.frame(dragonball = TRUE), jobname = "sponge_interaction_filter_singlesubset_2",
                    nodes = 1, cpus_per_node = 25)
results <- get_slurm_out(sjob, outtype = "raw")
results
