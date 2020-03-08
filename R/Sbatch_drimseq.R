library(rslurm)

drimseq <- function(dragonball){
  source("/nfs/home/students/evelyn/bachelor/R_workspace/MIRRETI/R/MIRRETI_drimseq.R")
  tpm.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/tcga_Kallisto_tpm"
  sampleannot.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/TCGA_phenotype_denseDataOnlyDownload.tsv"
  primary.disease <- "breast invasive carcinoma"
  ensembl_mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "http://jan2020.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")
  outfile <- paste("/nfs/home/students/evelyn/bachelor/data/thirdShot/datapacks/drimseq-1.0_datapack_",
                   gsub(" ", "", primary.disease), "_", 
                   gsub(" ", "", listEnsembl(ensembl_mart)$version[1]), ".RDS", sep = "")
  
  drimseq <- mirreti.drimseq(tpm.filepath = tpm.filepath, sampleannot.filepath = sampleannot.filepath, primary.disease = primary.disease, ensembl_mart = ensembl_mart, log = TRUE)
  saveRDS(drimseq, file = outfile)
}

sjob <- slurm_apply(drimseq, params = data.frame(dragonball = TRUE), jobname = "drimseq",
                    nodes = 1, cpus_per_node = 1)
