context('snpclust')
snpclust:::setup_temp_dir()

library(GWASTools)
gds <- snpclust:::save_hgdp_as_gds(check = TRUE)
bed <- snpclust:::gds_to_bed(gds)

ids <- 1:10
hgdp_gdata <- snpclust:::load_gds_as_genotype_data(gds)
test_haplo_mcmc <- function() {
  m_haplo <- snpclust:::haplo_mcmc(bed, getSnpVariable(hgdp_gdata, 'probe_id', ids))
  expect_is(m_haplo, 'matrix')
  expect_is(m_haplo[1], 'numeric')
  expect_true(!any(is.na(m_haplo)))
}
test_that('haplo_mcmc', test_haplo_mcmc())

suppressMessages(HGDP_GDATA_QC <- snpclust:::snprelate_qc(hgdp_gdata))

test_geno_block <- function() {
  geno = snpclust:::fetch_genotypes(hgdp_gdata, 1:10, 1:10, char = TRUE)
  expect_equal(dim(geno), c(10, 10))
}
test_that('geno_block', test_geno_block())

test_ld_select <- function() {
  ld = snpclust:::snprelate_ld_select(hgdp_gdata, snps_idx = 1:10,
    scans_idx = 1:10, min_r2 = 0.8) 
  expect_equal(length(ld$chr1), 9)
}
test_that('ld_select', test_ld_select())

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

snpclust_object <- snpclust::snpclust(gds = gds, n_axes = 3)
test_snpclust <- function() {
  expect_true(length(snpclust_object) == 8)
}
test_that('snpclust', test_snpclust())

context('ggplots')

test_ggplots <- function() {
  plts <- lapply(list(
      snpclust:::ggplot_pca(snpclust_object$pca, 'population',
        ellipses = TRUE),
      ggplot_manhat(pca = snpclust_object$pca, gdata = snpclust_object$gdata,
        axes = paste0('PC', 1:3)),
      ggplot_selection(peaks = snpclust_object$peaks,
        pca = snpclust_object$pca, axes = paste0('PC', 1:3))
    ), ggplot2::ggplotGrob) 
  lapply(plts, expect_is, 'gtable')
}
test_that('ggplots', test_ggplots())


