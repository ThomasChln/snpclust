
library(snpclust)
library(GWASTools)

paths <- paste0('../extdata/hgdp.', c('zip', 'txt', 'tar.gz', 'gds', 'rds'))
zippaths <- paste0('hgdp/', c('HGDP_FinalReport_Forward.txt', 'HGDP_Map.txt'))

.load_geno_data <- function(with_annotations = TRUE) {
  file <- system.file("extdata", "illumina_geno.nc", package = "GWASdata")
  nc <- GWASTools::NcdfGenotypeReader(file)

  snpAnnot <- scanAnnot <- NULL
  if (with_annotations) {
    data(illuminaSnpADF, package = "GWASdata")
    data(illuminaScanADF, package = "GWASdata")
    # need to rebuild old SNP annotation object to get allele methods
    snpAnnot <- SnpAnnotationDataFrame(pData(illuminaSnpADF))
    scanAnnot <- ScanAnnotationDataFrame(pData(illuminaScanADF))
  }

  GenotypeData(nc, snpAnnot = snpAnnot, scanAnnot = scanAnnot)
}

.get_test_gds_scan <- function(...) {
  genot <- load_gds_as_genotype_data(SNPRelate::snpgdsExampleFileName(), ...)

  nbscans <- nscan(genot)
  pheno <- data.frame(scanID = getScanID(genot),
    MYRESP = as.numeric(gl(2, k = 9, length = nbscans)) -1,
    MYCOV  = rnorm(nbscans, mean = 3, sd = 1), stringsAsFactors = FALSE
  )
  genot@scanAnnot <- ScanAnnotationDataFrame(pheno)
  validObject(genot)
  genot
}

GET_GDATA <- memoise::memoise(function() {
  .get_test_gds_scan(read_snp_annot = TRUE)
})

HGDP_DATA_PATH <- paths[4]

GET_HGDP_GDATA <- memoise::memoise(function() {
  snpclust::load_gds_as_genotype_data(HGDP_DATA_PATH)
})

GET_HGDP_GENO <- memoise::memoise(function() {
  t(snpclust::fetch_genotypes(GET_HGDP_GDATA()))
})

GET_HGDP_GDATA_QCED <- memoise::memoise(function() {
    snpclust::get_qced_data(GET_HGDP_GENO(), GET_HGDP_GDATA())
})

