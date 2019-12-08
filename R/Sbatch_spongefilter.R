library(rslurm)
library(SPONGE)


#---------------------------------------------------------------------------
#                               Main methode
#---------------------------------------------------------------------------


spongefilter <- function(preprocessed_data, coefficient.direction, out_file){
  load(as.character(preprocessed_data))
  mir_expr <- as.matrix(mir_expr)
  tpm_expr <- as.matrix(tpm_expr)
  interactions_matrix <- as.matrix(interactions_matrix)
  interaction_candidates_filtered <- sponge_gene_miRNA_interaction_filter(gene_expr = tpm_expr,
                                                                          mir_expr = mir_expr,
                                                                          mir_predicted_targets = interactions_matrix,
                                                                          coefficient.direction = as.character(coefficient.direction))
  saveRDS(interaction_candidates_filtered, file = as.character(out_file))
}



#---------------------------------------------------------------------------
#                               Parameters
#---------------------------------------------------------------------------


spongefilter.parameters <- as.data.frame(read.csv("/nfs/home/students/evelyn/bachelor/R_workspace/MIRRETI/params/SPONGEfilter_filetable.tsv", sep = "\t", header = T))


#---------------------------------------------------------------------------
#                             Slurm job request
#---------------------------------------------------------------------------


sjob <- slurm_apply(spongefilter, params = spongefilter.parameters, jobname = "sponge interaction filter",
                    nodes = 1, cpus_per_node = 40)
results <- get_slurm_out(sjob, outtype = "raw")
results
