library(rslurm)
library(SPONGE)
source("/nfs/home/students/evelyn/bachelor/R_workspace/MIRRETI/R/data_preprocess.R")


#---------------------------------------------------------------------------
#                               Main methode
#---------------------------------------------------------------------------

mirreti.runtime <- function(){
  mir.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena"
  tpm.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/tcga_Kallisto_tpm"
  utr5.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_5UTR.txt"
  cds.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_CDS.txt"
  utr3.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_3UTR.txt"
  sampleannot.filepath <- "/nfs/home/students/evelyn/bachelor/data/core_data/TCGA_phenotype_denseDataOnlyDownload.tsv"
  primary.disease <- "breast invasive carcinoma"
  sample.type.id = "11"
  
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
  
  
}