library(rslurm)


#---------------------------------------------------------------------------
#                               Main methode
#---------------------------------------------------------------------------


spongefilter <- function(dragonball){
  set.seed(01101990)
  source("/nfs/home/students/evelyn/bachelor/R_workspace/MIRRETI/R/MIRRETI.R")
  mir.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena"
  tpm.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/tcga_Kallisto_tpm"
  utr5.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_5UTR.txt"
  cds.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_CDS.txt"
  utr3.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_3UTR.txt"
  sampleannot.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/TCGA_phenotype_denseDataOnlyDownload.tsv"
  primary.disease <- "breast invasive carcinoma"
  sample.type.id = "01"
  ensembl_mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "http://apr2018.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")

  data <- mirreti(mir.filepath = mir.filepath, 
                  tpm.filepath = tpm.filepath,
                  interactions.filepath = list(utr5.filepath, cds.filepath, utr3.filepath),
                  sampleannot.filepath = sampleannot.filepath,
                  primary.disease = primary.disease,
                  sample.type.id = sample.type.id,
                  cluster.size = 25,
                  ensembl_mart = ensembl_mart)
  candidates <- data[[1]]
  interactions <- data[[2]]
  mir_expr <- data[[3]]
  tpm_expr <- data[[4]]
  
  save(candidates, interactions, mir_expr, tpm_expr, file = "/nfs/home/students/evelyn/bachelor/data/thirdShot/mirreti_datapack_breastcancer_primarytumor_vers92.RDS")
}



#---------------------------------------------------------------------------
#                             Slurm job request
#---------------------------------------------------------------------------


sjob <- slurm_apply(spongefilter, params = data.frame(dragonball = TRUE), jobname = "mirreti_ensemblversion92",
                    nodes = 1, cpus_per_node = 25)
results <- get_slurm_out(sjob, outtype = "raw")
results
