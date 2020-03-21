library(rslurm)
cluster.size = 10



#----------------------------------------------------------------------------------------------
#                                         Main method
#----------------------------------------------------------------------------------------------

sbatch.mirreti.4.1 <- function(primary.diseases, conditions, f.test, cluster.size){
   set.seed(01101990)
   log = TRUE
         start.time <- 0
         end.time <- 0
   
         if(log){
            cat("\n\t\t\t\t*** Welcome to MIRRETI 4.1 ***\n\n\n")
            cat("______________________________________________________________________\n\n")
            start.time <- Sys.time()
            cat(paste(start.time, "\tStart sourcing MIRRETI 4.1 methods and loading core data...\n", sep = ""))
         }
   source("/nfs/home/students/evelyn/bachelor/R_workspace/MIRRETI/R/MIRRETI_4.1.R")
   load("~/bachelor/data/core_data/mirreti_coredata.RDS")
   outfile.datapack <- paste("/nfs/home/students/evelyn/bachelor/data/datapacks/mirreti-4.1_datapack_", paste(gsub(" ", "", primary.diseases), collapse = "."), "_", paste(conditions, collapse = "."), "_Ens99.RDS", sep = "")
   outfile.minipack <- paste("/nfs/home/students/evelyn/bachelor/data/datapacks/mirreti-4.1_minipack_", paste(gsub(" ", "", primary.diseases), collapse = "."), "_", paste(conditions, collapse = "."), "_Ens99.RDS", sep = "")
         if(log){
            end.time <- Sys.time()
            cat(paste(end.time, "\t\U2705 \U1F1F8\U1F1FA\U1F1E8\U1F1E8\U1F1EA\U1F1F8\U1F1F8\U1F1EB\U1F1FA\U1F1F1\n", sep = ""))
            cat(paste("\U231B Time required:\t", round(difftime(end.time, start.time, units = "secs"), digits = 3), " secs\n\n", sep = ""))
         } 
   
   
   #-----------------------------------------------------------------------------------
   #                                Data preprocession
   #-----------------------------------------------------------------------------------
   
   conditions <- unlist(strsplit(as.character(conditions), "\\|"))
   primary.diseases <- unlist(strsplit(as.character(primary.diseases), "\\|"))
   data <- tcgaSamples.subset(tpm_expr = tpm_expr,
                              mir_expr = mir_expr,
                              sample_annotation = sample_annotation,
                              primary.diseases = primary.diseases,
                              sample.type.ids = conditions,
                              log = log)
   tpm_expr <- data[[1]]
   mir_expr <- data[[2]]
   sample_annotation <- data[[3]]
     rm(data)
     gc()
   interactions <- mirwalkInteractions.idconversion(interactions = interactions)
   
   data <- mirretiData.filter.actualTargets(tpm_expr = tpm_expr,
                                            mir_expr = mir_expr,
                                            interactions = interactions,
                                            var.threshold = 0.2)
   tpm_expr <- data[[1]]
   mir_expr <- data[[2]]
   interactions <- data[[3]]
     rm(data)
     gc()
   
     
   #-----------------------------------------------------------------------------------
   #                          One & Two condition MIRRETI
   #-----------------------------------------------------------------------------------
     
   # One condition
   if(length(primary.diseases) == 1 & length(conditions) == 1){
      interactions <- mirreti.oneCondition(tpm_expr = tpm_expr,
                                           mir_expr = mir_expr,
                                           interactions = interactions,
                                           f.test = f.test,
                                           cluster.size = cluster.size,
                                           log = log)
      save(tpm_expr, mir_expr, interactions, file = outfile.datapack)
   }
     
   # Two condition
   else{
      cat("START MIRRETI 4.1 TWO CONDITION ANALYSIS\n")
      interactions <- mirreti.twoConditions(tpm_expr = tpm_expr,
                                    mir_expr = mir_expr,
                                    sample_annotation = sample_annotation,
                                    interactions = interactions,
                                    cluster.size = cluster.size,
                                    log = log)
      #interactions <- data[[1]]
      #drimseq <- data[[2]]
      #   rm(data)
      #   gc()
      #save(tpm_expr, mir_expr, sample_annotation, interactions, drimseq, file = outfile.datapack)
      #save(interactions, drimseq, file = outfile.minipack)
      save(tpm_expr, mir_expr, sample_annotation, interactions, file = outfile.datapack)
   }
}



#----------------------------------------------------------------------------------------------
#                                      Server requests
#----------------------------------------------------------------------------------------------


# ONE CONDITION: breast cancer - primary tumor
params_BRCA_1 <- data.frame(primary.diseases = "breast invasive carcinoma", 
                            conditions = "1",
                            cluster.size = cluster.size)
sjob_BRCA_1 <- slurm_apply(sbatch.mirreti.4.1, 
                           params = params_BRCA_1, 
                           jobname = "mirreti_4.1_BRCA_1",
                           nodes = 1, 
                           cpus_per_node = cluster.size)
rm(params_BRCA_1, sjob_BRCA_1)
gc()



# TWO CONDITIONS: breast cancer - primary tumor & solid tissue
params_BRCA_1.11 <- data.frame(primary.diseases = "breast invasive carcinoma", 
                               conditions = "1|11",
                               cluster.size = cluster.size)
sjob_BRCA_1.11 <- slurm_apply(sbatch.mirreti.4.1, 
                              params = params_BRCA_1.11, 
                              jobname = "mirreti_4.1_BRCA_1.11",
                              nodes = 1, 
                              cpus_per_node = cluster.size)
rm(params_BRCA_1.11, sjob_BRCA_1.11)
gc()



# TWO CONDITIONS: breast cancer & kidney clear cell carcinoma - primary tumor
params_BRCA.KIRC_1 <- data.frame(primary.diseases = "breast invasive carcinoma|kidney clear cell carcinoma", 
                                 conditions = "1",
                                 cluster.size = cluster.size)
sjob_BRCA.KIRC_1 <- slurm_apply(sbatch.mirreti.4.1, 
                           params = params_BRCA.KIRC_1,
                           jobname = "mirreti_4.1_BRCA.KIRC_1",
                           nodes = 1, 
                           cpus_per_node = cluster.size)
rm(params_BRCA.KIRC_1, sjob_BRCA.KIRC_1)
gc()



# DIFFERENT CANCER TYPES ONE CONDITION ANALYSIS
params_KIRC_1 <- data.frame(primary.diseases = "kidney clear cell carcinoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_KIRC_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_KIRC_1, jobname = "mirreti_4.1_KIRC_1", nodes = 1, cpus_per_node = cluster.size)

params_LUAD_1 <- data.frame(primary.diseases = "lung adenocarcinoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_LUAD_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_LUAD_1, jobname = "mirreti_4.1_LUAD_1", nodes = 1, cpus_per_node = cluster.size)

params_OV_1 <- data.frame(primary.diseases = "ovarian serous cystadenocarcinoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_OV_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_OV_1, jobname = "mirreti_4.1_OV_1", nodes = 1, cpus_per_node = cluster.size)

params_LUSC_1 <- data.frame(primary.diseases = "lung squamous cell carcinoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_LUSC_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_LUSC_1, jobname = "mirreti_4.1_LUSC_1", nodes = 1, cpus_per_node = cluster.size)

params_GBM_1 <- data.frame(primary.diseases = "glioblastoma multiforme", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_GBM_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_GBM_1, jobname = "mirreti_4.1_GBM_1", nodes = 1, cpus_per_node = cluster.size)

params_HNSC_1 <- data.frame(primary.diseases = "head & neck squamous cell carcinoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_HNSC_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_HNSC_1, jobname = "mirreti_4.1_HNSC_1", nodes = 1, cpus_per_node = cluster.size)

params_UCEC_1 <- data.frame(primary.diseases = "uterine corpus endometrioid carcinoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_UCEC_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_UCEC_1, jobname = "mirreti_4.1_UCEC_1", nodes = 1, cpus_per_node = cluster.size)

params_STAD_1 <- data.frame(primary.diseases = "stomach adenocarcinoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_STAD_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_STAD_1, jobname = "mirreti_4.1_STAD_1", nodes = 1, cpus_per_node = cluster.size)

params_THCA_1 <- data.frame(primary.diseases = "thyroid carcinoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_THCA_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_THCA_1, jobname = "mirreti_4.1_THCA_1", nodes = 1, cpus_per_node = cluster.size)

params_PRAD_1 <- data.frame(primary.diseases = "prostate adenocarcinoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_PRAD_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_PRAD_1, jobname = "mirreti_4.1_PRAD_1", nodes = 1, cpus_per_node = cluster.size)

params_COAD_1 <- data.frame(primary.diseases = "colon adenocarcinoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_COAD_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_COAD_1, jobname = "mirreti_4.1_COAD_1", nodes = 1, cpus_per_node = cluster.size)

params_LGG_1 <- data.frame(primary.diseases = "brain lower grade glioma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_LGG_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_LGG_1, jobname = "mirreti_4.1_LGG_1", nodes = 1, cpus_per_node = cluster.size)

params_SKCM_1 <- data.frame(primary.diseases = "skin cutaneous melanoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_SKCM_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_SKCM_1, jobname = "mirreti_4.1_SKCM_1", nodes = 1, cpus_per_node = cluster.size)

params_LIHC_1 <- data.frame(primary.diseases = "liver hepatocellular carcinoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_LIHC_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_LIHC_1, jobname = "mirreti_4.1_LIHC_1", nodes = 1, cpus_per_node = cluster.size)

params_BLCA_1 <- data.frame(primary.diseases = "bladder urothelial carcinoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_BLCA_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_BLCA_1, jobname = "mirreti_4.1_BLCA_1", nodes = 1, cpus_per_node = cluster.size)

params_KIRP_1 <- data.frame(primary.diseases = "kidney papillary cell carcinoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_KIRP_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_KIRP_1, jobname = "mirreti_4.1_KIRP_1", nodes = 1, cpus_per_node = cluster.size)

params_CESC_1 <- data.frame(primary.diseases = "cervical & endocervical cancer", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_CESC_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_CESC_1, jobname = "mirreti_4.1_CESC_1", nodes = 1, cpus_per_node = cluster.size)

params_SARC_1 <- data.frame(primary.diseases = "sarcoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_SARC_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_SARC_1, jobname = "mirreti_4.1_SARC_1", nodes = 1, cpus_per_node = cluster.size)

params_ESCA_1 <- data.frame(primary.diseases = "esophageal carcinoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_ESCA_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_ESCA_1, jobname = "mirreti_4.1_ESCA_1", nodes = 1, cpus_per_node = cluster.size)

params_LAML_1 <- data.frame(primary.diseases = "acute myeloid leukemia", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_LAML_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_LAML_1, jobname = "mirreti_4.1_LAML_1", nodes = 1, cpus_per_node = cluster.size)

params_PAAD_1 <- data.frame(primary.diseases = "pancreatic adenocarcinoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_PAAD_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_PAAD_1, jobname = "mirreti_4.1_PAAD_1", nodes = 1, cpus_per_node = cluster.size)

params_PCPG_1 <- data.frame(primary.diseases = "pheochromocytoma & paraganglioma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_PCPG_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_PCPG_1, jobname = "mirreti_4.1_PCPG_1", nodes = 1, cpus_per_node = cluster.size)

params_READ_1 <- data.frame(primary.diseases = "rectum adenocarcinoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_READ_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_READ_1, jobname = "mirreti_4.1_READ_1", nodes = 1, cpus_per_node = cluster.size)

params_TGCT_1 <- data.frame(primary.diseases = "testicular germ cell tumor", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_TGCT_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_TGCT_1, jobname = "mirreti_4.1_TGCT_1", nodes = 1, cpus_per_node = cluster.size)

params_THYM_1 <- data.frame(primary.diseases = "thymoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_THYM_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_THYM_1, jobname = "mirreti_4.1_THYM_1", nodes = 1, cpus_per_node = cluster.size)

params_ACC_1 <- data.frame(primary.diseases = "adrenocortical cancer", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_ACC_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_ACC_1, jobname = "mirreti_4.1_ACC_1", nodes = 1, cpus_per_node = cluster.size)

params_KICH_1 <- data.frame(primary.diseases = "kidney chromophobe", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_KICH_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_KICH_1, jobname = "mirreti_4.1_KICH_1", nodes = 1, cpus_per_node = cluster.size)

params_MESO_1 <- data.frame(primary.diseases = "mesothelioma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_MESO_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_MESO_1, jobname = "mirreti_4.1_MESO_1", nodes = 1, cpus_per_node = cluster.size)

params_UVM_1 <- data.frame(primary.diseases = "uveal melanoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_UVM_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_UVM_1, jobname = "mirreti_4.1_UVM_1", nodes = 1, cpus_per_node = cluster.size)

params_UCS_1 <- data.frame(primary.diseases = "uterine carcinosarcoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_UCS_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_UCS_1, jobname = "mirreti_4.1_UCS_1", nodes = 1, cpus_per_node = cluster.size)

params_DLBC_1 <- data.frame(primary.diseases = "diffuse large B-cell lymphoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_DLBC_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_DLBC_1, jobname = "mirreti_4.1_DLBC_1", nodes = 1, cpus_per_node = cluster.size)

params_CHOL_1 <- data.frame(primary.diseases = "cholangiocarcinoma", conditions = "1", f.test = "both", cluster.size = cluster.size)
sjob_CHOL_1 <- slurm_apply(sbatch.mirreti.4.1, params = params_CHOL_1, jobname = "mirreti_4.1_CHOL_1", nodes = 1, cpus_per_node = cluster.size)

rm(params_CHOL_1, sjob_CHOL_1, params_DLBC_1, sjob_DLBC_1, params_UCS_1, sjob_UCS_1, params_UVM_1, sjob_UVM_1, params_MESO_1, sjob_MESO_1, params_KICH_1, sjob_KICH_1, params_ACC_1, sjob_ACC_1, params_THYM_1, sjob_THYM_1, params_TGCT_1, sjob_TGCT_1, params_READ_1, sjob_READ_1, params_PCPG_1, sjob_PCPG_1, params_PAAD_1, sjob_PAAD_1, params_LAML_1, sjob_LAML_1, params_ESCA_1, sjob_ESCA_1, params_SARC_1, sjob_SARC_1, params_CESC_1, sjob_CESC_1, params_KIRP_1, sjob_KIRP_1, params_BLCA_1, sjob_BLCA_1, params_LIHC_1, sjob_LIHC_1, params_SKCM_1, sjob_SKCM_1, params_LGG_1, sjob_LGG_1, params_COAD_1, sjob_COAD_1, params_PRAD_1, sjob_PRAD_1, params_THCA_1, sjob_THCA_1, params_STAD_1, sjob_STAD_1, params_UCEC_1, sjob_UCEC_1, params_HNSC_1, sjob_HNSC_1, params_GBM_1, sjob_GBM_1, params_LUSC_1, sjob_LUSC_1, params_OV_1, sjob_OV_1, params_LUAD_1, sjob_LUAD_1, params_KIRC_1, sjob_KIRC_1)
gc()