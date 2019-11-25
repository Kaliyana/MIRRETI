library(rslurm)


#---------------------------------------------------------------------------
#                               Main methode
#---------------------------------------------------------------------------

exprformater <- function(expr_file, sample_annot_file, primary_disease, sample_type_id, out_file){
  source("/nfs/home/students/evelyn/bachelor/R_workspace/MIRRETI/R/Expression_file_formater.R")
  
  expr_data <- readExpressionFile(as.character(expr_file))
  sample_annotation <- read.csv(as.character(sample_annot_file), header = T, sep = "\t")
  expr_data_subset <- process_expression_file(expr_data = expr_data,
                                              sample_annotation = sample_annotation,
                                              primary_disease = as.character(primary_disease),
                                              sample_type_id = as.character(sample_type_id))
  saveRDS(expr_data_subset, file = as.character(out_file))
}


#---------------------------------------------------------------------------
#                               Parameters
#---------------------------------------------------------------------------

exprformater.parameters <- as.data.frame(read.csv(paste(getwd(), "params", "Exprformater_filetable.tsv", sep = .Platform$file.sep), header = T, sep = "\t"))


#---------------------------------------------------------------------------
#                             Slurm job request
#---------------------------------------------------------------------------

sjob <- slurm_apply(exprformater, params = exprformater.parameters, jobname = "exprfile formating",
                    nodes = 1, cpus_per_node = 4)
results <- get_slurm_out(sjob, outtype = "raw")
results
