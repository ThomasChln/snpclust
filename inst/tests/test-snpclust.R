
library(snpclust)
library(GWASTools)

gds <- snpclust:::save_hgdp_as_gds()
tar <- snpclust:::gds_to_bedtargz(gds)

context('Dim. red.')

ids <- 1:10
hgdp_gdata <- snpclust:::load_gds_as_genotype_data(gds)
test_haplo_mcmc <- function() {
  m_haplo <- snpclust::haplo_mcmc(tar, getSnpVariable(hgdp_gdata, 'probe_id', ids))
  expect_is(m_haplo, 'matrix')
  expect_is(m_haplo[1], 'numeric')
  expect_true(!any(is.na(m_haplo)))
}
test_that('haplo_mcmc', test_haplo_mcmc())

context('qc pca')

suppressMessages(HGDP_GDATA_QC <- snpclust:::snprelate_qc(hgdp_gdata))

test_snprelate_qc <- function() {
  expect_equal(length(HGDP_GDATA_QC), 2)
  expect_identical(names(HGDP_GDATA_QC), c('gdata', 'df_qc'))
  expect_is(HGDP_GDATA_QC$gdata, 'GenotypeData')
}
test_that('snprelate_qc', test_snprelate_qc())

test_snprelate_pca <- function() {
  test <- snpclust:::snprelate_pca(HGDP_GDATA_QC$gdata)
  expect_is(test, c("data.frame", "qb_pcafort"))
}
test_that('snprelate_pca', test_snprelate_pca())
close(hgdp_gdata)

context('snpclust')

snpclust_object <- snpclust::snpclust(gds = gds)
test_snpclust <- function() {
  expect_true(length(snpclust_object) == 6)
  file.remove(gds)
}
test_that('snpclust', test_snpclust())

context('ggplots')

test_ggplots <- function() {
  plts <- lapply(list(
      ggplot_pca(snpclust_object$pca, 'population', ellipses = TRUE),
      ggplot_manhat(pca = snpclust_object$pca, gdata = snpclust_object$gdata),
      ggplot_selection(peaks = snpclust_object$peaks, pca = snpclust_object$pca)
    ), ggplot2::ggplotGrob) 
  lapply(append(plts, plot_pca_pairs(1:3, snpclust_object$pca, 'population')),
    expect_is, 'gtable')
}
test_that('ggplots', test_ggplots())

context('misc qc')

test_sample_impute <- function() {
  num_mat <- matrix(1:16, 4)
  num_mat[2:3, ] <- matrix(rep(NA, 8), 2)
  num_mat <- sample_impute(num_mat)
  apply(num_mat, 2, function(column) {
      expect_true(all(column[2:3] %in% column[c(1,4)]))
    })
  char_mat <- matrix(letters[1:16], 4)
  char_mat[2:3, ] <- matrix(rep(NA, 8), 2)
  char_mat <- sample_impute(char_mat)
  apply(char_mat, 2, function(column) {
      expect_true(all(column[2:3] %in% column[c(1,4)]))
    })
}
test_that('sample_impute', test_sample_impute())

test_merge_dfs <- function() {
  dfs <- lapply(1:3, function(i) iris[sample(1:150, 100), ])
  expected <- merge(merge(dfs[[1]], dfs[[2]]), dfs[[3]])
  expect_identical(expected, merge_dfs(dfs))
}
test_that('merge_dfs', test_merge_dfs())

.df_rbind_all <- function() {
  one <- mtcars[1:4, ]
  two <- mtcars[11:14, ]

  expected <- mtcars[c(1:4, 11:14), ]
  rownames(expected) <- NULL

  expect_identical(df_rbind_all(one, two), expected)
  expect_is(df_rbind_all(one, two, one), "data.frame")
  expect_identical(df_rbind_all(list(one, two)), expected)


  expect_identical(df_rbind_all(one, two, use_row_names = TRUE),
     mtcars[c(1:4, 11:14), ])

  expect_identical(df_rbind_all(one, two, use_row_names = TRUE),
     mtcars[c(1:4, 11:14), ])

  ### factors
  df <- df_rbind_all(iris, iris)
  expect_true(is.factor(df$Species))

  df1 <- data.frame(x = as.factor(letters[1:2]))
  df <- df_rbind_all(df1, df1)
  expect_true(is.factor(df$x))

  ### different levels
  df2 <- data.frame(x = as.factor(letters[3:4]))
  df <- df_rbind_all(df1, df2)
  expect_false(is.factor(df$x))

  ### included levels
  df3 <- data.frame(x = as.factor(letters[2]))
  df <- df_rbind_all(df1, df3)
  expect_false(is.factor(df$x))
}
test_that('df_rbind_all', .df_rbind_all())


context('geno_data')

test_actg_tsv_to_gdata  <- function() {
  dir <- paste0(gsub('tests$', '', getwd()), 'extdata')
  paths <- paste0(dir, '/hgdp.', c('zip', 'txt'))
  zippaths <- paste0('hgdp/', c('HGDP_FinalReport_Forward.txt', 'HGDP_Map.txt'))
  txts_paths <- unzip(paths[1], zippaths, exdir = dir, junkpaths = TRUE)
  actg_gdata <- snpclust:::actg_tsv_to_gdata(txts_paths[1], paths[2],
    c('scan_id', 'gender', 'population', 'geographic_origin', 'region'),
    txts_paths[2])

  expect_is(actg_gdata, 'GenotypeData')

  # Open all
  file_paths <- c(paths[2], txts_paths[2])
  l_files <- list(data.table::fread(txts_paths[1], '\t'))
  l_files[2:3] <- lapply(file_paths, data.table::fread, '\t')
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

test_bed_targz_to_gds <- function() {
  gdspath <- SNPRelate::snpgdsExampleFileName()
  snpclust:::setup_temp_dir()
  names <- paste0('example_data', c('.gds', '_test.gds'))

  # tidy example file for plink: sex and pheno as 2 first slots of sample annot
  file.copy(gdspath, names[1])
  snpclust:::.format_gds_for_BED2GDS(
    names[1],
    'sample.annot',
    'sex',
    'pop.group'
  )

  # gds -> gdata
  gdata <- snpclust:::load_gds_as_genotype_data(names[1])
  on.exit(close(gdata))

  # gdata -> bed targz
  snpclust:::save_genotype_data_as_plink(gdata, names[2], verbose = FALSE)
  tar(names[2], list.files(), compression = 'gzip', tar = Sys.getenv('TAR'))

  # bed targz -> gdata 2
  snpclust:::bed_targz_to_gds(names[2], names[2])
  gdata_test <- snpclust:::load_gds_as_genotype_data(names[2])
  on.exit(close(gdata_test))

  # gdata == gdata 2

  # geno and snp annot
  getters <- list(snpclust:::fetch_genotypes, getSnpID, getChromosome, getPosition,
    getAlleleA, getAlleleB)
  lapply(getters, function(getter) {
      l_data_and_test <- lapply(list(gdata, gdata_test), getter)
      expect_identical(l_data_and_test[[1]], l_data_and_test[[2]])
    })

  # sample annot
  plink_cols <- paste0('sample.annot/', c('sex', 'phenotype'))
  l_gdatas <- list(gdata, gdata_test)
  l_plink_infos <- lapply(l_gdatas, function(gdata) {
      gds <- snpclust:::request_snpgds_file(gdata)$snpgds
      l_plink_info <- lapply(plink_cols, snpclust:::.gds_read_index, gds)
    })
  sapply(1:2, function(info_idx) {
      expect_identical(l_plink_infos[[1]][[info_idx]],
        l_plink_infos[[1]][[info_idx]])
    })
}
test_that('bed_targz_to_gds', test_bed_targz_to_gds())


