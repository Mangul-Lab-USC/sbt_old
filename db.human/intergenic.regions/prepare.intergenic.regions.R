library(TxDb.Hsapiens.UCSC.hg19.knownGene)
genic <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
genic <- reduce(genic, ignore.strand=T)
intergenic <- gaps(genic)
intergenic <- intergenic[strand(intergenic) == "*"]
write.table(genes_df, "intergenic.regions.hg19.txt", append = FALSE, sep = " ", dec = ".",row.names = TRUE, col.names = TRUE)
