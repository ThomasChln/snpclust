
source('helper_setup.R')
context("GenotypeDataSubset")

my_openfn.gds <- function(filename) {
  gdsfmt::openfn.gds(filename, TRUE, TRUE, TRUE)
}


.load_sample_geno_data <- function() {
  data(affyScanADF, package = 'GWASdata')
  data(affySnpADF, package = 'GWASdata')
  file <- system.file("extdata", "affy_geno.nc", package="GWASdata")
  nc <- GWASTools::NcdfGenotypeReader(file)
  gdata <- GenotypeData(nc, snpAnnot = affySnpADF, scanAnnot = affyScanADF)

  set.seed(666)
  snpdf <- Biobase::pData(gdata@snpAnnot)
  snpdf$alleleA <- sample(c('C', 'A', 'T', '-'), nrow(snpdf), replace = TRUE)
  snpdf$alleleB <- sample(c('-', 'T', 'A', 'G'), nrow(snpdf), replace = TRUE)
#  snpdf$alleleA[1] <- NA
#  snpdf$alleleB[1] <- NA
  Biobase::pData(gdata@snpAnnot) <- snpdf

  gdata
}


.gdata_are_identical <- function(g2, gdata) {
  is(g2, 'GenotypeData') && is (gdata, 'GenotypeData') &&
  identical(getSnpID(g2), getSnpID(gdata)) &&
    identical(getScanID(g2), getScanID(gdata)) &&
    identical(getGenotype(g2), getGenotype(gdata))
}

.test_sub_gdata <- function(g, gdata, snps, scans) {
  expect_is(g, 'GenotypeDataSubset')
  expect_equal(nsnp(g), length(snps))
  expect_equal(nscan(g), length(scans))
  expect_equal(getSnpID(g), getSnpID(gdata, snps))
  expect_equal(getSnpID(g, 1), getSnpID(gdata, snps[1]))

  expect_equal(getScanID(g), getScanID(gdata, scans))
  expect_equal(getScanID(g, 1), getScanID(gdata, scans[1]))

  expect_equal(hasSnpAnnotation(g),
    hasSnpAnnotation(gdata))
  expect_equal(hasScanAnnotation(g),
     hasScanAnnotation(gdata))
  expect_equal(getChromosome(g), getChromosome(gdata, snps))
  expect_equal(getPosition(g), getPosition(gdata, snps))

  expect_equal(getAlleleA(g), getAlleleA(gdata, snps))
  expect_equal(getAlleleB(g), getAlleleB(gdata, snps))
  expect_equal(hasSex(g), hasSex(gdata))

  suppressWarnings(expect_equal(getSex(g), getSex(gdata)))

  suppressWarnings(expect_equal(getSnpVariableNames(g),
      getSnpVariableNames(gdata)))

  suppressWarnings(vars <- getSnpVariableNames(gdata))
  for (v in vars) {
    expect_equal(hasSnpVariable(g, v),
      hasSnpVariable(gdata, v))
    expect_equal(getSnpVariable(g, v)
    , getSnpVariable(gdata, v))
  }

  suppressWarnings(expect_equal(getScanVariableNames(g),
      getScanVariableNames(gdata)))
  suppressWarnings(vars <- getScanVariableNames(gdata))
  for (v in vars) {
    expect_equal(hasScanVariable(g, v),
      hasScanVariable(gdata, v))
    suppressWarnings(expect_equal(getScanVariable(g, v),
        getScanVariable(gdata, v)))
  }

  expect_identical(getGenotype(g), getGenotype(gdata)[snps, scans])
  expect_identical(getGenotype(g, snp = c(1, 1), scan = c(1,1)),
    getGenotype(gdata)[snps[1], scans[1]])


}

################################### TESTS ####################################

.get_scan_annot <- function() {
  snps <- 10:99
  scans <- 35:43

  GDATA <- GET_GDATA()

  g2 <- genotype_data_subset(GDATA, snps, scans)
  sa <- get_scan_annot(g2)
  expect_is(sa, 'ScanAnnotationDataFrame')
  expect_equivalent(nrow(sa), length(scans))
  expect_identical(getScanID(sa), getScanID(g2))
  expect_identical(getScanID(sa), getScanID(GDATA)[scans])
}
test_that("get_scan_annot", .get_scan_annot())



.get_snp_annot <- function() {
  snps <- 10:99
  scans <- 35:43

  gdata <- .load_geno_data(TRUE)

  g2 <- genotype_data_subset(gdata, snps, scans)
  sa <- get_snp_annot(g2)

  expect_is(sa, 'SnpAnnotationDataFrame')
  expect_equivalent(nrow(sa), length(snps))
  expect_identical(getSnpID(sa), getSnpID(g2))
  expect_identical(getSnpID(sa), getSnpID(gdata)[snps])
}
test_that("get_snp_annot", .get_snp_annot())



.io <- function() {
  snps <- 10:99
  scans <- 35:43
  gdata <- .load_geno_data(TRUE)
  g1 <- genotype_data_subset(gdata, snps, scans)

  setup_temp_dir()

  save_genotype_data_as_gds(g1, 'g1.snpgds')
  g2 <- load_gds_as_genotype_data('g1.snpgds')

  expect_equivalent(class(g1), "GenotypeDataSubset")
  expect_equivalent(class(g2), "GenotypeData")

  expect_true(.gdata_are_identical(g2, g1))

  expect_that(check_genotype_data(g2), not(throws_error()))

}
test_that("io", .io())

.genotype_data_subset <- function() {
  gdata <- .load_geno_data(TRUE)

  gdata0 <- gdata
  gdata0@snpAnnot <- NULL
  gdata0@scanAnnot <- NULL

  snps <- 1000L
  scans <- 10:50

  g2 <- genotype_data_subset(gdata, snps, scans)
  .test_sub_gdata(g2, gdata, snps, scans)
  g3 <- genotype_data_subset(gdata0, snps, scans)
  .test_sub_gdata(g3, gdata0, snps, scans)

  snps <- 10:99
  scans <- 1L
  g2 <- genotype_data_subset(gdata, snps, scans)
  .test_sub_gdata(g2, gdata, snps, scans)
  g3 <- genotype_data_subset(gdata0, snps, scans)
  .test_sub_gdata(g3, gdata0, snps, scans)

  snps <- 10:99
  scans <- 35:43
  g2 <- genotype_data_subset(gdata, snps, scans)
  .test_sub_gdata(g2, gdata, snps, scans)
  g3 <- genotype_data_subset(gdata0, snps, scans)
  .test_sub_gdata(g3, gdata0, snps, scans)

  ## check blocks
  snps <- seq.int(10L, 99L, by = 3L)
  scans <- seq.int(35L, 50L, by = 2L)
  g2 <- genotype_data_subset(gdata, snps, scans)

  m1 <- getGenotype(g2, snp = c(3L, 5L), scan = c(4L, 3L))
  m2 <- fetch_genotypes(gdata, snps[seq.int(3, length.out = 5)],
     scans[seq.int(4, length.out = 3)])
  expect_identical(m1, m2)

  m1 <- getGenotype(g2, snp = c(7, -1), scan = c(4, -1))
  m2 <- fetch_genotypes(gdata, snps[7:length(snps)], scans[4:length(scans)])
  expect_identical(m1, m2)

  snps <- 12L
  scans <- 7L
  g2 <- genotype_data_subset(gdata, snps, scans)
  .test_sub_gdata(g2, gdata, snps, scans)
  g3 <- genotype_data_subset(gdata0, snps, scans)
  .test_sub_gdata(g3, gdata0, snps, scans)

  g2 <- genotype_data_subset(gdata)
  expect_true(.gdata_are_identical(g2, gdata))

  ### subset of subset
  snps <- 10:13
  scans <- 35:43
  g2 <- genotype_data_subset(gdata, snps, scans)
  g3 <- genotype_data_subset(g2, 1:2, 1:3)
  expect_identical(getSnpID(g3, 1), getSnpID(gdata, snps[1]))
  expect_identical(getSnpID(g3),  getSnpID(gdata, snps[1:2]))
  expect_identical(getScanID(g3), getScanID(gdata, scans[1:3]))

}
test_that("test_genotype_data_subset", .genotype_data_subset())

context("fetch_genotypes")

GDATA <- GET_GDATA()

.fetch_genotypes_by_chunk_iterator <- function() {
  nb_snps <- nsnp(GDATA)

  .test_it <- function(it) {
    total <- 0L
    for (i in seq_len(it$nb)) {
      chunk <- it$get(i)
      expect_identical(chunk, fetch_genotypes(GDATA,
          snps_idx = seq.int(from = total + 1, length.out = nrow(chunk))))
      total <- total + nrow(chunk)
    }
    expect_equal(total, nb_snps)
  }

  it <- fetch_genotypes_by_chunk_iterator(GDATA, chunk_size = 997)

  .test_it(it)
  expect_equal(it$chunk_size, 997)

  it <- fetch_genotypes_by_chunk_iterator(GDATA, nb_chunks = 13)
  expect_equal(it$nb, 13)
  .test_it(it)

  expect_error(it$get(0), 'bad chunk index')
  expect_error(it$get(14), 'bad chunk index')

}
test_that("fetch_genotypes_by_chunk_iterator_genotypes",
  .fetch_genotypes_by_chunk_iterator())




.fetch_genotypes <- function() {
  ### default params
  genos <- getGenotype(GDATA)
  expect_identical(fetch_genotypes(GDATA), genos)

  ### edge cases
  expect_error(fetch_genotypes(GDATA, 0:5), 'invalid')
  expect_error(fetch_genotypes(GDATA, nsnp(GDATA) + 1), 'invalid')
  expect_error(fetch_genotypes(GDATA, scans_idx = 0:5), 'invalid')
  expect_error(fetch_genotypes(GDATA, scans_idx = nscan(GDATA) + 1),
    'invalid')

  ### snps_idx
  # one snp
  expect_identical(fetch_genotypes(GDATA, 1), genos[1, , drop = FALSE])
  expect_identical(fetch_genotypes(GDATA, nsnp(GDATA)),
    genos[nsnp(GDATA), , drop = FALSE])

  # interval
  expect_identical(fetch_genotypes(GDATA, 5:10), genos[5:10, ])
  # order
  expect_identical(fetch_genotypes(GDATA, 22:17), genos[22:17, ])

  ### both
  idx1 <- c(seq.int(10, 20, 2), seq.int(21, 11, -2), 50:40)
  idx2 <- c(2:8, seq.int(10, 20, 2), seq.int(41, 31, -2))
  expect_identical(fetch_genotypes(GDATA, idx1, idx2, snps_first = TRUE),
    genos[idx1, idx2])
  expect_identical(fetch_genotypes(GDATA, idx1, idx2, snps_first = FALSE),
    genos[idx1, idx2])

  # with subset
  g2 <- genotype_data_subset(GDATA, 1L , 1L)
  genos <- fetch_genotypes(g2)
  expect_identical(genos, fetch_genotypes(GDATA, 1, 1))

}
test_that("fetch_genotypes", .fetch_genotypes())




.fetch_genotypes_by_idx <- function() {
  gdata <- GDATA
  genos <- getGenotype(gdata)

  # by block
  idx <- c(seq.int(10, 20, 2), seq.int(21, 11, -2))
  # 2 blocks, but restore order
  expect_identical(fetch_genotypes_by_idx(gdata, idx, snps_first = TRUE),
     genos[idx, ])
  expect_identical(fetch_genotypes_by_idx(gdata, idx, snps_first = FALSE),
     genos[idx, ])

   expect_identical(fetch_genotypes_by_idx(gdata, scans_idx = idx),
     genos[, idx])


  ### snps_idx
  # one snp
  expect_identical(fetch_genotypes_by_idx(gdata, 1), genos[1, , drop = FALSE])
  expect_identical(fetch_genotypes_by_idx(gdata, nsnp(gdata)),
    genos[nsnp(gdata), , drop = FALSE])

  # interval
  expect_identical(fetch_genotypes_by_idx(gdata, 5:10), genos[5:10, ])
  # order
  expect_identical(fetch_genotypes_by_idx(gdata, 22:17), genos[22:17, ])

  ### scans_idx
  # one scan
  expect_identical(fetch_genotypes_by_idx(gdata, , 1), genos[, 1, drop = FALSE])
  expect_identical(fetch_genotypes_by_idx(gdata, , nscan(gdata),
       snps_first = TRUE), genos[, nscan(gdata), drop = FALSE])
  expect_identical(fetch_genotypes_by_idx(gdata, , nscan(gdata),
     snps_first = FALSE), genos[, nscan(gdata), drop = FALSE])

  # interval
  expect_identical(fetch_genotypes_by_idx(gdata, , 5:10), genos[, 5:10])
  # order
  expect_identical(fetch_genotypes_by_idx(gdata, , 22:17), genos[, 22:17])

  # by block
  idx <- c(seq.int(10, 20, 2), seq.int(21, 11, -2))
  # 1 block, but restore order
  expect_identical(fetch_genotypes_by_idx(gdata, , idx), genos[, idx])

  ### both
  idx1 <- c(seq.int(10, 20, 2), seq.int(21, 11, -2), 50:40)
  idx2 <- c(2:8, seq.int(10, 20, 2), seq.int(41, 31, -2))
  expect_identical(fetch_genotypes_by_idx(gdata, idx1, idx2, snps_first = TRUE),
     genos[idx1, idx2])
  expect_identical(fetch_genotypes_by_idx(gdata, idx1, idx2, snps_first = FALSE),
     genos[idx1, idx2])
}
test_that("fetch_genotypes_by_idx", .fetch_genotypes_by_idx())



.fetch_genotypes_by_blocks <- function() {
  gdata <- GDATA

  ### single blocks
  snps_blocks = list(c(5, 9))
  scans_blocks = list(c(21, 17))
  ref <- getGenotype(gdata, snp = snps_blocks[[1]], scan = scans_blocks[[1]])

  genos <- fetch_genotypes_by_blocks(gdata, snps_blocks, scans_blocks,
    snps_first = TRUE)
	expect_identical(genos, ref)
  genos <- fetch_genotypes_by_blocks(gdata, snps_blocks, scans_blocks,
    snps_first = FALSE)
	expect_identical(genos, ref)

  ### multiple blocks
  snps_blocks = list(c(100, 99), c(5, 9), c(1000, 5))
  scans_blocks = list(c(21, 17), c(3, 4), c(100, 9))
  m11 <- getGenotype(gdata, snp = snps_blocks[[1]], scan = scans_blocks[[1]])
  m12 <- getGenotype(gdata, snp = snps_blocks[[1]], scan = scans_blocks[[2]])
  m13 <- getGenotype(gdata, snp = snps_blocks[[1]], scan = scans_blocks[[3]])
  m21 <- getGenotype(gdata, snp = snps_blocks[[2]], scan = scans_blocks[[1]])
  m22 <- getGenotype(gdata, snp = snps_blocks[[2]], scan = scans_blocks[[2]])
  m23 <- getGenotype(gdata, snp = snps_blocks[[2]], scan = scans_blocks[[3]])
  m31 <- getGenotype(gdata, snp = snps_blocks[[3]], scan = scans_blocks[[1]])
  m32 <- getGenotype(gdata, snp = snps_blocks[[3]], scan = scans_blocks[[2]])
  m33 <- getGenotype(gdata, snp = snps_blocks[[3]], scan = scans_blocks[[3]])
  m1 <- cbind(m11, m12, m13)
  m2 <- cbind(m21, m22, m23)
  m3 <- cbind(m31, m32, m33)
  ref <- rbind(m1, m2, m3)
  genos <- fetch_genotypes_by_blocks(gdata, snps_blocks, scans_blocks,
    snps_first = TRUE)
  expect_identical(genos, ref)
  genos <- fetch_genotypes_by_blocks(gdata, snps_blocks, scans_blocks,
    snps_first = FALSE)
  expect_identical(genos, ref)

  ### singletons (check drop)
  snps_blocks = list(c(5, 1))
  scans_blocks = list(c(21, 1))
  ref <- getGenotype(gdata, snp = snps_blocks[[1]], scan = scans_blocks[[1]])
  genos <- fetch_genotypes_by_blocks(gdata, snps_blocks, scans_blocks,
    snps_first = TRUE)
  expect_equal(dim(genos), c(1L, 1L))
  expect_identical(genos[1, 1], ref)

  genos2 <- fetch_genotypes_by_blocks(gdata, snps_blocks, scans_blocks,
    snps_first = FALSE)
  expect_equal(genos2, genos)
}
test_that("fetch_genotypes_by_blocks", .fetch_genotypes_by_blocks())

context('input/output')

.save_genotype_data_as_plink <- function() {
  dir <- qbutils::setup_temp_dir(chdir = FALSE)
  #  GDS based
  save_genotype_data_as_plink(GDATA, 'test1', dir, snps_idx = 20:30,
      scans_idx = 5:10, verbose = FALSE)

  paths <- file.path(dir, paste0('test1', c('.bed', '.fam', '.bim')))
  expect_true(all(file.exists(paths)))
  gdsfn <- file.path(dir, 'test1.gds')
  snpgdsBED2GDS(paths[1], paths[2], paths[3], gdsfn, verbose = FALSE)

  gdata <- load_gds_as_genotype_data(gdsfn)
  expect_identical(fetch_genotypes(gdata), fetch_genotypes(GDATA, 20:30, 5:10))
  expect_identical(fetch_genotypes(gdata, char = TRUE),
      fetch_genotypes(GDATA, 20:30, 5:10, char = TRUE))

  # no GDS already there
  g2 <- genotype_data_subset(GDATA, snps_idx = 20:30, scans_idx = 5:10)
  save_genotype_data_as_plink(g2, 'test2', dir, verbose = FALSE)
  paths <- file.path(dir, paste0('test2', c('.bed', '.fam', '.bim')))
  gdsfn <- file.path(dir, 'test2.gds')
  snpgdsBED2GDS(paths[1], paths[2], paths[3], gdsfn, verbose = FALSE)
  g3 <- load_gds_as_genotype_data(gdsfn)
  expect_identical(fetch_genotypes(g3, char = TRUE),
      fetch_genotypes(g2, char = TRUE))
}
test_that("save_genotype_data_as_plink", .save_genotype_data_as_plink())

.io_gds_as_genotype_data <- function() {
  getAlleleA <- getAlleleA
  getAlleleB <- getAlleleB


  gdata <- .load_sample_geno_data()
  check_genotype_data(gdata)
  output <- tempfile()
  on.exit({unlink(output)})

  save_genotype_data_as_gds(gdata, output)

  # check alleles
  gds <- my_openfn.gds(output)

  alleles <- sort(unique(snpclust:::read_snp_gds_alleles(gds)))
  all_ref <- paste(getAlleleA(gdata), getAlleleB(gdata), sep = '/')
  all_ref[is.na(getAlleleA(gdata)) | is.na(getAlleleB(gdata))] <- NA
  all_ref <- sort(unique(all_ref))
  expect_identical(alleles, all_ref)
  gdsfmt::closefn.gds(gds)

  gdata2 <- load_gds_as_genotype_data(output)
  check_genotype_data(gdata2)

  expect_identical(getGenotype(gdata2, char = TRUE),
    getGenotype(gdata, char = TRUE))
  expect_identical(getGenotype(gdata2), getGenotype(gdata))

  expect_identical(gdata2@snpAnnot, gdata@snpAnnot)
  expect_identical(gdata2@scanAnnot, gdata@scanAnnot)
  close(gdata2)

  gdata3 <- load_gds_as_genotype_data(output, read_scan_annot = FALSE,
    read_snp_annot = FALSE)
  check_genotype_data(gdata3)
  expect_identical(getGenotype(gdata3), getGenotype(gdata))
  expect_identical(gdata3@snpAnnot, NULL)
  expect_identical(gdata3@scanAnnot, NULL)
  close(gdata3)

  save_genotype_data_as_gds(gdata, output, save_annotations = FALSE,
    check = FALSE)
  gdata2 <- load_gds_as_genotype_data(output)

  expect_identical(fetch_genotypes(gdata2), fetch_genotypes(gdata))
  expect_identical(fetch_genotypes(gdata2, char = TRUE),
    fetch_genotypes(gdata, char = TRUE))
#  expect_identical(getGenotype(gdata2, char = TRUE),
#    getGenotype(gdata, char = TRUE))

}
test_that("load and save gds_as_genotype_data", .io_gds_as_genotype_data())


