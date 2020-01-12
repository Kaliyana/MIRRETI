source("/nfs/home/students/evelyn/bachelor/R_workspace/MIRRETI/R/MIRRETI.R")
library(profvis)


#-------------------------------------------------------------------------------------
#                                   Core data
#-------------------------------------------------------------------------------------


# rows: mirs           cols: samples
# 743 obs. 10824 variables
mir_expr <- read.expression.file("/nfs/home/students/evelyn/bachelor/data/core_data/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena")

# rows: transcripts    cols: samples
# 197044 obs. 10535 variables
tpm_expr <- read.expression.file("/nfs/home/students/evelyn/bachelor/data/core_data/tcga_Kallisto_tpm")
  
# $miRNA: mirs         $mRNA: transcripts      $Genesymbol: genesymbol
# 110466891 obs. 6 variables
interactions <- read.interaction.files("/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_5UTR.txt",
                                       "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_CDS.txt",
                                       "/nfs/home/students/evelyn/bachelor/data/core_data/hsa_miRWalk_3UTR.txt")
  
# sample    sample_type_id    sample_type   X_primary_disease
# 12804 obs. 4 variables
sample_annotation <- read.sample.annotation.file("/nfs/home/students/evelyn/bachelor/data/core_data/TCGA_phenotype_denseDataOnlyDownload.tsv")



#-------------------------------------------------------------------------------------
#                                   Preprocession
#-------------------------------------------------------------------------------------


mirreti.breastinvasivecarcinoma <- function(mir_expr, tpm_expr, interactions, sample_annotation){
  interactions <- preproc.idconversion.refseqtoensembl(interactions)  # 116.529.799 obs. 7 variables
  
  box <- preproc.filter.samples(mir_expr, tpm_expr, sample_annotation, primary.disease = "breast invasive carcinoma", sample.type.id = "11")
  mir_expr <- box[[1]]                                # 743 obs. 90 variables
  tpm_expr <- box[[2]]                                # 197.044 obs. 113 variables
  
  mir_expr <- preproc.filter.na(mir_expr, 0)          # 485 obs. 90 variables
  tpm_expr <- preproc.filter.na(tpm_expr, -9.9658)    # 84.103 obs. 112 variables
  
  mir_expr <- preproc.filter.variance(mir_expr)       # 388 obs. 90 variables
  tpm_expr <- preproc.filter.variance(tpm_expr)       # 67.282 obs. 112 variables
  
  box <- preproc.crossfilter(mir_expr, tpm_expr, interactions)
  mir_expr <- box[[1]]                                # 378 obs. 89 variables
  tpm_expr <- box[[2]]                                # 15.834 obs. 89 variables
  interactions <- box[[3]]                            # 5.027.525 obs. 7 variables
  
  box <- preproc.data.spongefilter(mir_expr, tpm_expr, interactions)
  mir_expr <- box[[1]]                                # num [1:89, 1:378]
  tpm_expr <- box[[2]]                                # num [1:89, 1:15.834]
  interactions_matrix <- box[[3]]                     # 15.834 obs. 378 variables
}

