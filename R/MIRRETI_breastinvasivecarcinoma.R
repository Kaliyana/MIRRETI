source("/nfs/home/students/evelyn/bachelor/R_workspace/MIRRETI/R/data_preprocess.R")

#rows: mirs           cols: samples
mir_expr <- read.expression.file("/nfs/home/students/evelyn/bachelor/data/TCGA_expression_data/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena")
#rows: transcripts    cols: samples
tpm_expr <- read.expression.file("/nfs/home/students/evelyn/bachelor/data/TCGA_expression_data/tcga_Kallisto_tpm")
#$miRNA: mirs         $mRNA: transcripts      $Genesymbol: genesymbol
utr3 <- read.interaction.file("/nfs/home/students/evelyn/bachelor/data/hsa_miRWalk_interactions/hsa_miRWalk_3UTR.txt")
utr3$genomic_region <- rep("3UTR", nrow(utr3))
utr5 <- read.interaction.file("/nfs/home/students/evelyn/bachelor/data/hsa_miRWalk_interactions/hsa_miRWalk_5UTR.txt")
utr5$genomic_region <- rep("5UTR", nrow(utr5))
cds <- read.interaction.file("/nfs/home/students/evelyn/bachelor/data/hsa_miRWalk_interactions/hsa_miRWalk_CDS.txt")
cds$genomic_region <- rep("CDS", nrow(cds))
#interactions <- read.interaction.file("/nfs/home/students/evelyn/bachelor/data/hsa_miRWalk_interactions/hsa_miRWalk_wholeGene.txt")
interactions <- rbind(utr5, cds, utr3)

sample_annotation <- read.sample.annotation.file("/nfs/home/students/evelyn/bachelor/data/TCGA_expression_data/cancertype_filter/TCGA_phenotype_denseDataOnlyDownload.tsv")

outfile <- "/nfs/home/students/evelyn/bachelor/data/core_data/mirreti_coredata.RDS"
saveRDS(mir_expr, tpm_expr, interactions, sample_annotation, file = outfile)