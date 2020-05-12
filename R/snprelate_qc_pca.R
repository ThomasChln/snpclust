
#' snprelate_qc
#'
#' Quality control using SNPRelate functions
#'
#' @param gdata       Genotype data object
#' @param samples_nas NA threshold for samples, default 3 pct
#' @param ibs         Samples identity by state threshold, default 99 pct
#' @param keep_ids    Samples ids to keep even if IBS is higher than threshold.
#'                    Used for monozygotic twins.
#' @param snps_nas    NA threshold for SNPs, default 1 pct
#' @param maf         Minor allele frequency threshold, default 5 pct
#' @param tagsnp      TagSNP r2 correlation threshold, default 0.8
#' @param n_cores     Number of cores
#' @return List of gdata, Genotype data object, and df_qc, QC info data frame
snprelate_qc <- function(gdata, samples_nas = 0.03, ibs = 0.99, keep_ids = NULL,
  snps_nas = 0.01, maf = 0.05, tagsnp = 0.8, n_cores = 1) {

  stopifnot(inherits(gdata, 'GenotypeData'))
  samples = getScanID(gdata)
  snps = getSnpID(gdata)
  df_qc <- .rbind_qc(NULL, 'Raw', NA, samples, snps)
  gds <- request_snpgds_file(gdata)$snpgds

  # samples NAs
  nas_rates <- suppressMessages(
    SNPRelate::snpgdsSampMissRate(gds, samples, snps))
  samples <- samples[nas_rates <= samples_nas]
  df_qc <- .rbind_qc(df_qc, 'Samples NAs', samples_nas, samples, snps)

  # samples IBS
  m_ibs <- suppressMessages(SNPRelate::snpgdsIBS(gds, samples, snps,
      num.thread = n_cores, verbose = FALSE))$ibs
  m_ibs[lower.tri(m_ibs, TRUE)] <- NA
  max_ibs <- apply(m_ibs[, -1], 2, max, na.rm = TRUE)
  samples_ibs <- samples[c(TRUE, max_ibs < ibs)]

  # keep ids
  if (!is.null(keep_ids) &&
    any(keep_ids %in% samples)) {
    samples_ibs <- union(samples_ibs, keep_ids[keep_ids %in% samples])
  }
  samples <- samples_ibs
  df_qc <- .rbind_qc(df_qc, 'Identity by state - Twins', ibs, samples, snps)

  # SNPs NAs
  snps <- unlist(suppressMessages(SNPRelate::snpgdsSelectSNP(gds,
        samples, snps, missing.rate = snps_nas, verbose = FALSE)))
  df_qc <- .rbind_qc(df_qc, 'SNPs NAs', snps_nas, samples, snps)

  # MAF
  snps <- unlist(suppressMessages(SNPRelate::snpgdsSelectSNP(gds,
       samples, snps, maf = maf, verbose = FALSE)))
  df_qc <- .rbind_qc(df_qc, 'MAF', maf, samples, snps)

  # TagSNP
  if (!is.na(tagsnp)) {
    snps <- unlist(suppressMessages(SNPRelate::snpgdsLDpruning(gds,
          samples, snps, ld.threshold = sqrt(tagsnp), method = 'r',
          verbose = FALSE)))
    df_qc <- .rbind_qc(df_qc, 'TagSNP', tagsnp, samples, snps)
  }
  samples = match(samples, getScanID(gdata))
  snps = match(snps, getSnpID(gdata))
  gdata <- genotype_data_subset(gdata, snps, samples)

  list(gdata = gdata, qc = df_qc)
}

#' snprelate_pca
#'
#' Principal component analysis using SNPRelate
#'
#' @param gdata   Genotype data object
#' @param n_axes  Number of axes
#' @param n_cores Number of cores
#' @return Data frame with variance, observations, and variables
snprelate_pca <- function(gdata, n_axes = 32, n_cores = 1) {
  stopifnot(inherits(gdata, 'GenotypeData'))

  # get pca
  gds <- request_snpgds_file(gdata)$snpgds
  samples = getScanID(gdata)
  snps = getSnpID(gdata)
  n_axes <- min(c(nscan(gdata), nsnp(gdata)) - 1, n_axes)
  pca <- suppressMessages(SNPRelate::snpgdsPCA(gds, samples, snps,
    num.thread = n_cores, eigen.cnt = n_axes, verbose = FALSE))
  if (n_axes == 0) n_axes <- ncol(pca$eigenvect)
  seq_axes <- seq_len(n_axes)
  name_axes <- paste0('PC', seq_axes)

  # get variance
  var <- pca$eigenval[seq_axes] / sum(pca$eigenval, na.rm = TRUE)
  var_percent <- round(var * 100)
  df_variance <- data.frame(matrix(c(var, var_percent), 2, byrow = TRUE),
    DIMRED_VARNAME = c('Explained_variance', 'Explained_variance_percent'),
    DIMRED_VARTYPE = 'OTHER')
  names(df_variance)[seq_axes] <- name_axes

  # get obs
  mat_obs <- pca$eigenvect * sqrt(length(samples))
  df_obs <- data.frame(mat_obs, DIMRED_VARTYPE = 'OBS',
    DIMRED_VARNAME = paste0('OBS_', samples), stringsAsFactors = FALSE)
  names(df_obs)[seq_axes] <- name_axes

  # get vars
  df_vars <- data.frame(DIMRED_VARNAME = paste0('VAR_', snps),
    DIMRED_VARTYPE = 'VAR', stringsAsFactors = FALSE)
  m_vars <- suppressMessages(SNPRelate::snpgdsPCASNPLoading(pca, gds, n_cores,
      FALSE))
  df_vars <- cbind(t(m_vars$snploading), df_vars)
  names(df_vars)[seq_axes] <- name_axes

  df_pca <- df_rbind_all(df_variance, df_obs, df_vars)

  # add obs annotations from gdata
  obs_idxs <- df_pca$DIMRED_VARTYPE == 'OBS'
  df_scans_annots = gdata_scans_annots(gdata)
  df_annot <- subset(df_scans_annots, scanID %in% samples) 
  for (colname in names(df_scans_annots)) {
    df_pca[[colname]] <- NA
    df_pca[obs_idxs, colname] <- as.character(df_annot[[colname]])
    df_pca[[colname]] <- factor(df_pca[[colname]])
  }
  
  class(df_pca) <- c('data.frame', 'pca')

  df_pca
}

.rbind_qc <- function(df_qc, Step, Parameter, samples, snps) {
  df_qc_new_line <- data.frame(Step, Parameter,
    Samples = length(samples), SNPs = length(snps))
  rbind(df_qc, df_qc_new_line)
}

sample_impute <- function(data) {
  stopifnot(is.matrix(data))
  na_idxs <- which(colSums(is.na(data)) != 0)
  data[, na_idxs] <- apply(data[, na_idxs, drop = FALSE], 2,
    function(variable) {
      nas <- is.na(variable)
      variable[nas] <- sample(variable[!nas], sum(nas), TRUE)
      variable
    })

  data
}


merge_dfs <- function(l_df, ...) {
  stopifnot(is.list(l_df))
  while (length(l_df) > 1) {
    l_df[[2]] <- merge(l_df[[1]], l_df[[2]], ...)
    l_df <- l_df[-1]
  }

  l_df[[1]]
}


df_rbind_all <- function(...,  use_row_names = FALSE) {
  df <- as.data.frame(suppressWarnings(dplyr::bind_rows(...)))
  if (use_row_names) {
    dots <- list(...)
    dfs <-  if (is.list(dots[[1]]) && !is.data.frame(dots[[1]]))
      dots[[1]] else dots
    row_names <- unlist(lapply(dfs, rownames))
    rownames(df) <- row_names
  }

  df
}

snprelate_ld_select <- function(gdata,
  window_length = 500L,
  min_r2,
  window_size = NA,
  snps_idx = NULL,
  scans_idx = NULL,
  remove.monosnp = FALSE,
  autosome.only = FALSE,
  method = 'r',
  threads = 1,
  quiet = FALSE,
  ...
) {
  info <- if (quiet) function(...) {} else function(...) message(sprintf(...))

  info('managing input')
  gdsinfo <- request_snpgds_file(gdata, snps_idx, scans_idx)
  gds <- gdsinfo$snpgds
  if (gdsinfo$new_file) {
    on.exit({closefn.gds(gds); unlink(gds$filename)}, TRUE)
  }

  info('SNPRelate::snpgdsLDpruning')

  # there is currently a bug in snpgdsLDpruning with huge values
  bp <- min(1e9, window_length * 1000)

  res <- SNPRelate::snpgdsLDpruning(gds, gdsinfo$scan_ids, gdsinfo$snp_ids,
    slide.max.bp = bp,
    slide.max.n = window_size,
    ld.threshold = sqrt(min_r2),
    remove.monosnp = remove.monosnp,
    autosome.only = autosome.only,
    method = method, num.thread = threads,
    verbose = !quiet)

  res
}

snps_hla <- function(gdata) {
  location <- getSnpVariable(gdata, c('chromosome', 'position'))
  in_hla <- location$chromosome == 6
  in_hla <- in_hla & location$position > 25e6 & location$position < 35e6
}

# Character vector for subsetting the dataset:
#                'hla' use only SNPs in HLA
#                'non_hla' use only SNPs outside HLA
#                'no_controls' exclude samples with Control in phenotype 
#                If several subsets are provided, each slot will be a list
#                with slots corresponding to each subsets.
gdata_subsets <- function(subsets, gdata) {
  snps_idx <- switch(set,
    hla = which(snps_hla(gdata)),
    non_hla = which(!snps_hla(gdata)),
    GWASTools::getSnpID(gdata))

  df_scans_annots = gdata_scans_annots(gdata)
  if (set == 'no_controls' && !'phenotype' %in% names(df_scans_annots)) {
    set <- ''
  }
  scans_idx <- switch(set,
    no_controls = which(df_scans_annots$phenotype != 'Control'),
    seq_len(GWASTools::nscan(gdata)))

  genotype_data_subset(gdata, snps_idx, scans_idx)
}
