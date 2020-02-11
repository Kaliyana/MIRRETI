library(rslurm)


#---------------------------------------------------------------------------
#                               Main methode
#---------------------------------------------------------------------------


mirreti <- function(dragonball){
  set.seed(01101990)
  source("/nfs/home/students/evelyn/bachelor/R_workspace/MIRRETI/R/MIRRETI.R")
  mir.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena"
  tpm.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/tcga_Kallisto_tpm"
  utr5.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_5UTR.txt"
  cds.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_CDS.txt"
  utr3.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_3UTR.txt"
  sampleannot.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/TCGA_phenotype_denseDataOnlyDownload.tsv"
  primary.disease <- "breast invasive carcinoma"
  sample.type.id = "11"
  
  #ensembl_mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "http://apr2018.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")     # vers 92
  ensembl_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")                                                               # vers 99
  
  outfile <- paste("/nfs/home/students/evelyn/bachelor/data/thirdShot/candidate/events_with_tpms/mirreti_datapack_",
                   gsub(" ", "", primary.disease), "_", 
                   sample.type.id, "_",
                   gsub(" ", "", listEnsembl(ensembl_mart)$version[1]), ".RDS", sep = "")
  

  data <- mirreti(mir.filepath = mir.filepath, 
                  tpm.filepath = tpm.filepath,
                  interactions.filepath = list(utr5.filepath, cds.filepath, utr3.filepath),
                  sampleannot.filepath = sampleannot.filepath,
                  primary.disease = primary.disease,
                  sample.type.id = sample.type.id,
                  cluster.size = 25,
                  ensembl_mart = ensembl_mart)
  candidate_events <- data[[1]]
  tpm_expr <- data[[2]]
  
  save(candidate_events, tpm_expr, file = outfile)
}



#---------------------------------------------------------------------------
#                             Slurm job request
#---------------------------------------------------------------------------


sjob <- slurm_apply(mirreti, params = data.frame(dragonball = TRUE), jobname = "mirreti_breastcancer_11_vers99",
                    nodes = 1, cpus_per_node = 25)
results <- get_slurm_out(sjob, outtype = "raw")
results
