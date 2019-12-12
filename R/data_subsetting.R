library(SPONGE)
library(Rfast)


filename <- "D:\\Bioinformatics\\Bachelordata\\secondShot\\subset_analysis\\testdataTop10_g_spongefiltered"

mir_expr <- mir_expr
tpm_expr <- as.matrix(subset.exprdata.for.variance(as.data.frame(gene_expr), 10))
interactions_matrix <- as.matrix(targetscan_symbol[rownames(targetscan_symbol) %in% colnames(tpm_expr), colnames(targetscan_symbol) %in% colnames(mir_expr)])


load("D:\\Bioinformatics\\Bachelordata\\secondShot\\spongerun_data_breastinvasivecarcinome_solidtissuenormal.RDS")
mir_expr <- as.matrix(mir_expr)
tpm_expr <- as.matrix(subset.exprdata.for.variance(tpm_expr, 100))
interactions_matrix <- as.matrix(interactions_matrix[rownames(interactions_matrix) %in% colnames(tpm_expr), ])

sponge.filter(tpm_expr = tpm_expr,
              mir_expr = mir_expr,
              interactions_matrix = interactions_matrix,
              outfile = filename,
              print = T)



sponge.filter <- function(tpm_expr, mir_expr, interactions_matrix, outfile, print = F){
  sponge_filtered_le <- sponge_gene_miRNA_interaction_filter(gene_expr = tpm_expr,
                                                             mir_expr = mir_expr,
                                                             mir_predicted_targets = interactions_matrix,
                                                             coefficient.direction = "<")
  sponge_filtered_ge <- sponge_gene_miRNA_interaction_filter(gene_expr = tpm_expr,
                                                             mir_expr = mir_expr,
                                                             mir_predicted_targets = interactions_matrix,
                                                             coefficient.direction = ">")
  
  save(sponge_filtered_le, sponge_filtered_ge, file = paste(filename, ".RDS", sep = ""))
  
  if(print){
    sink(paste(filename, "_le.txt", sep = ""))
    print(sponge_filtered_le)
    sink()
    
    sink(paste(filename, "_ge.txt", sep = ""))
    print(sponge_filtered_ge)
    sink()
  }
}
