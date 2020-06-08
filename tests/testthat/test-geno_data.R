context('genotype data')
snpclust:::setup_temp_dir()

test_actg_tsv_to_gdata  <- function() {
  paths <- system.file('extdata', paste0('hgdp.', c('zip', 'txt')),
    package = 'snpclust')
  zippaths <- file.path('hgdp',
    c('HGDP_FinalReport_Forward.txt', 'HGDP_Map.txt'))
  txts_paths <- unzip(paths[1], zippaths, junkpaths = TRUE)
  actg_gdata <- snpclust:::actg_tsv_to_gdata(txts_paths[1], paths[2],
    c('scan_id', 'gender', 'population', 'geographic_origin', 'region'),
    txts_paths[2])

  expect_is(actg_gdata, 'GenotypeData')

  # Open all
  file_paths <- c(paths[2], txts_paths[2])
  l_files <- list(data.table::fread(txts_paths[1], sep = '\t'))
  l_files[2:3] <- lapply(file_paths, data.table::fread, sep = '\t')
  l_files <- lapply(l_files, as.data.frame)
  file.remove(txts_paths)
  geno_test <- l_files[[1]][-1]
  scans_test <- l_files[[2]]
  snps_test <- l_files[[3]]

  # Check SNPs order
  snps <- actg_gdata@snpAnnot@data$probe_id
  # Order following chromosome and position
  snps_test[2] <- as.integer(factor(snps_test[[2]],
      c(1:22, 'X', 'XY', 'Y', 'M'), nmax = 26, exclude = NULL))
  snps_order <- order(snps_test[[2]], snps_test[[3]])
  snps_test <- snps_test[[1]][snps_order]
  expect_true(all(snps == snps_test))

  # Check scans order
  scans <- actg_gdata@snpAnnot@data$scan_id
  # Discard missings
  scans_test <- scans_test
  scans_test_merge <- match(scans_test[[1]], colnames(geno_test))
  scans_test <- scans_test[[1]][!is.na(scans_test_merge)]
  expect_true(all(scans == scans_test))

  geno <- snpclust:::fetch_genotypes(actg_gdata)
  geno[is.na(geno)] <- 3L
  geno_test <- geno_test[snps_order, ]
  geno_test <- geno_test[, na.omit(scans_test_merge)]
  geno_test <- t(apply(geno_test, 1, snpclust:::actg_to_numeric, '--'))
  geno_test[is.na(geno_test)] <- 3L
  expect_true(all(geno == geno_test))
}
test_that('Convert hgdp text files to gds', test_actg_tsv_to_gdata())


gdspath <- SNPRelate::snpgdsExampleFileName()
gdata <- snpclust:::load_gds_as_genotype_data(gdspath)

test_bed_to_gds <- function() {
  names <- paste0('example_data', c('.gds', '_test.gds'))

  # tidy example file for plink: sex and pheno as 2 first slots of sample annot
  file.copy(gdspath, names[1])
  snpclust:::.format_gds_for_BED2GDS(names[1], 'sample.annot', 'sex',
    'pop.group')

  # gdata -> bed 
  suppressWarnings(snpclust:::save_genotype_data_as_plink(gdata, 'plink_test',
    verbose = FALSE))

  # bed  -> gdata 2
  snpclust:::bed_to_gds(paste0('plink_test.', c('bed', 'bim', 'fam')),
    names[2])
  gdata_test <- snpclust:::load_gds_as_genotype_data(names[2])
  on.exit(close(gdata_test))

  # gdata == gdata 2
  snpclust:::check_genotype_data(gdata_test)

  # geno and snp annot
  getters <- list(snpclust:::fetch_genotypes, getSnpID, getChromosome, getPosition,
    getAlleleA, getAlleleB)
  lapply(getters, function(getter) {
      l_data_and_test <- lapply(list(gdata, gdata_test), getter)
      expect_identical(l_data_and_test[[1]], l_data_and_test[[2]])
    })

  # sample annot
  plink_cols <- 'sample.annot/sex'
  l_gdatas <- list(gdata, gdata_test)
  l_plink_infos <- lapply(l_gdatas, function(gdata) {
      gds <- snpclust:::request_snpgds_file(gdata)$snpgds
      l_plink_info <- lapply(plink_cols, snpclust:::.gds_read_index, gds)
    })
  expect_identical(l_plink_infos[[1]][[1]],
    l_plink_infos[[1]][[1]])
}
test_that('bed_to_gds', test_bed_to_gds())

close(gdata)

gds = snpclust:::save_hgdp_as_gds()
gdata = snpclust:::load_gds_as_genotype_data(gds)
on.exit(close(gdata))

test_plink_merge <- function() {
  library(magrittr)
  gdata1 = snpclust:::genotype_data_subset(gdata, 1:10, 1:10)
  snpclust:::save_genotype_data_as_gds(gdata1, 'tmp1.gds')
  gdata2 = snpclust:::genotype_data_subset(gdata, 6:15, 1:10)
  snpclust:::save_genotype_data_as_gds(gdata2, 'tmp2.gds')
  snpclust:::snprelate_gds_to_bed('tmp1')
  snpclust:::snprelate_gds_to_bed('tmp2')

  lapply(c('tmp1', 'tmp2'), paste0, '.', c('bed', 'bim', 'fam')) %>%
    snpclust:::plink_merge(., '.', 'plink_merge') %>%
    sapply(file.exists) %>% all %>% expect_true
}
test_that('plink_merge', test_plink_merge())
