library(rslurm)
library(SPONGE)


#---------------------------------------------------------------------------
#                               Main methode
#---------------------------------------------------------------------------


spongefilter <- function(mir_data, tpm_data, interaction_file, coefficient.direction, out_file){
  load(as.character(mir_data))
  load(as.character(tpm_data))
  interaction_candidates <- as.matrix(read.csv(as.character(interaction_file), sep = "\t", header = T, row.names = 1, check.names = F))
  tpm_expr <- as.matrix(extended_tpm_expr)
  colnames(extended_mirna_expr) <- gsub("\\.", "-", colnames(extended_mirna_expr))
  mirna_expr <- as.matrix(extended_mirna_expr)
  interaction_candidates_filtered <- sponge_gene_miRNA_interaction_filter(gene_expr = tpm_expr,
                                                                          mir_expr = mirna_expr,
                                                                          mir_predicted_targets = interaction_candidates,
                                                                          coefficient.direction = as.character(coefficient.direction))
  saveRDS(interaction_candidates_filtered, file = as.character(out_file))
}



#---------------------------------------------------------------------------
#                               Parameters
#---------------------------------------------------------------------------


spongefilter.parameters <- as.data.frame(read.csv("/nfs/home/students/evelyn/bachelor/R_workspace/OwnRpackage/R/SPONGEfilter_filetable.tsv", sep = "\t", header = T))



#---------------------------------------------------------------------------
#                             Slurm job request
#---------------------------------------------------------------------------


sjob <- slurm_apply(spongefilter, params = spongefilter.parameters, jobname = "sponge interaction filter",
                    nodes = 1, cpus_per_node = 40)
results <- get_slurm_out(sjob, outtype = "raw")
results
