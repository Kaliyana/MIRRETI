library(biomaRt)


load("D:\\Bioinformatics\\Bachelordata\\secondShot\\subset_analysis\\tcgadataTop10_t_spongefiltered.RDS")
ge_candidates <- unlist.sponge.candidates(sponge_filtered_ge)
le_candidates <- unlist.sponge.candidates(sponge_filtered_le)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)
filters <- listFilters(ensembl)

values <- unique(unlist(ge_candidates$transcript))
candidates_annot <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_exon_id', 'start_position', 'end_position', 'transcript_start', 'transcript_end', 'exon_chrom_start', 'exon_chrom_end'),
                          filters = 'ensembl_transcript_id',
                          values = values,
                          mart = ensembl)

unlist.sponge.candidates <- function(sponge_filtered){
  candidates <- do.call(rbind, sponge_filtered)
  candidates$transcript <- lapply(strsplit(rownames(candidates), "\\."), '[', 1)
  rownames(candidates) <- NULL
  return(candidates)
}
