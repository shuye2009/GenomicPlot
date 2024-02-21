# create data for examples and vignettes

gtffile <- system.file("extdata", "gencode.v19.annotation_chr19.gtf",
package = "GenomicPlot")

txdb <- custom_TxDb_from_GTF(gtffile, genome = "hg19")

AnnotationDbi::saveDb(txdb,
    "C:/GREENBLATT/Rscripts/GenomicPlot/inst/extdata/txdb.sql")

gf5_meta <- GenomicPlot::prepare_5parts_genomic_features(txdb, meta = TRUE,
    nbins = 100, fiveP = -2000, threeP = 1000, longest = TRUE)
usethis::use_data(gf5_meta, overwrite = TRUE)

gf5_genomic <- GenomicPlot::prepare_5parts_genomic_features(txdb, meta = FALSE,
    nbins = 100, fiveP = -2000, threeP = 1000, longest = TRUE)
usethis::use_data(gf5_genomic, overwrite = TRUE)
