library(biomaRt)
library(dplyr)


#-------------------------------------------------------------------------------------
#               biomaRt ID conversion RefSeq mRNA -> Ensembl transcripts
#-------------------------------------------------------------------------------------

get_ensemblids <- function(refseqids){
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  refseq_ensembl <- getBM(attributes = c('refseq_mrna', 'ensembl_transcript_id', 'external_gene_name'),
                          filters = 'refseq_mrna',
                          values = refseqids,
                          mart = ensembl)
  predrefseq_ensembl <- getBM(attributes = c('refseq_mrna_predicted', 'ensembl_transcript_id', 'external_gene_name'),
                              filters = 'refseq_mrna_predicted',
                              values = refseqids,
                              mart = ensembl)
  colnames(predrefseq_ensembl)[1] <- "refseq_mrna"
  combined_mappings <- rbind(refseq_ensembl, predrefseq_ensembl)
  return(combined_mappings)
}


dispose_blank_mappings <- function(refseq_ensembl){
  refseq_ensembl <- refseq_ensembl[!(refseq_ensembl$ensembl_transcript_id == ""),]
  return(refseq_ensembl)
}

extract_distinct_mappings <- function(refseq_ensembl){
  double_entries <- unique(refseq_ensembl$refseq_mrna[duplicated(refseq_ensembl$refseq_mrna)])
  distinct_mappings <- refseq_ensembl[!(refseq_ensembl$refseq_mrna %in% double_entries),]
  return(distinct_mappings)
}

extract_ambiguous_mappings <- function(refseq_ensembl){
  double_entries <- unique(refseq_ensembl$refseq_mrna[duplicated(refseq_ensembl$refseq_mrna)])
  ambiguous_mappings <- refseq_ensembl[refseq_ensembl$refseq_mrna %in% double_entries,]
  return(ambiguous_mappings)
}




#-------------------------------------------------------------------------------------
#                                     EXAMPLE RUN
#-------------------------------------------------------------------------------------

build_interaction_matrix <- function(interaction_candidates, refseq_ensembl_mappings){

}


#-------------------------------------------------------------------------------------
#                                     EXAMPLE RUN
#-------------------------------------------------------------------------------------

interaction_candidates <- read.csv("D:\\Bioinformatics\\Bachelordata\\hsa_miRWalk_interactions\\hsa_miRWalk_5UTR.txt", sep = "\t", header = T)
refseq_ensembl_mapping <- get_ensemblids(unique(interaction_candidates$mRNA))
refseq_ensembl_noblanks <- dispose_blank_mappings(refseq_ensembl_mapping)
refseq_ensembl_distincts <- extract_distinct_mappings(refseq_ensembl_mapping)
refseq_ensembl_ambiguous <- extract_ambiguous_mappings(refseq_ensembl_mapping)


