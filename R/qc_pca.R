
.rbind_qc <- function(df_qc, Step, Parameter, l_ids) {
  df_qc_new_line <- data.frame(Step, Parameter, lapply(l_ids, length))
  names(df_qc_new_line)[3:4] <- c('Samples', 'SNPs')
  rbind(df_qc, df_qc_new_line)
}

snprelate_qc <- function(gdata, sample_nas = .03, snp_nas = .01, maf = .05,
  tsnp = .8, ibs = .99, monozygotic_twins_ids = NULL, n_cores = 2) {

  stopifnot(inherits(gdata, 'GenotypeData'))
  l_ids_original <- list(getScanID(gdata), getSnpID(gdata))
  l_ids <- l_ids_original
  df_qc <- .rbind_qc(NULL, 'Raw', NA, l_ids)
  gds <- request_snpgds_file(gdata)$snpgds

  # samples NAs
  sample_nas_rates <- SNPRelate::snpgdsSampMissRate(gds, l_ids[[1]], l_ids[[2]])
  l_ids[[1]] <- l_ids[[1]][sample_nas_rates <= sample_nas]
  df_qc <- .rbind_qc(df_qc, 'Samples NAs', sample_nas, l_ids)

  # samples IBS
  m_ibs <- SNPRelate::snpgdsIBS(gds, l_ids[[1]], l_ids[[2]], num.thread = n_cores,
    verbose = FALSE)$ibs
  m_ibs[lower.tri(m_ibs, TRUE)] <- NA
  max_ibs <- apply(m_ibs[, -1], 2, max, na.rm = TRUE)
  l_ids[[1]] <- l_ids[[1]][c(TRUE, max_ibs < ibs)]

  # monozygotic twins
  if (!is.null(monozygotic_twins_ids) &&
    any(monozygotic_twins_ids %in% l_ids_original[[1]])) {
    l_ids[[1]] <- union(l_ids[[1]],
      monozygotic_twins_ids[monozygotic_twins_ids %in% l_ids_original[[1]]])
  }
  df_qc <- .rbind_qc(df_qc, 'Identity by state - Twins', ibs, l_ids)

  # MAF
  l_ids[[2]] <- unlist(SNPRelate::snpgdsSelectSNP(gds, l_ids[[1]], l_ids[[2]],
      missing.rate = snp_nas, verbose = FALSE))
  df_qc <- .rbind_qc(df_qc, 'SNPs NAs', snp_nas, l_ids)
  l_ids[[2]] <- unlist(SNPRelate::snpgdsSelectSNP(gds, l_ids[[1]], l_ids[[2]],
      maf = maf, verbose = FALSE))
  df_qc <- .rbind_qc(df_qc, 'MAF', maf, l_ids)

  # TagSNP
  if (!is.na(tsnp)) {
    l_ids[[2]] <- unlist(SNPRelate::snpgdsLDpruning(gds, l_ids[[1]], l_ids[[2]],
        ld.threshold = sqrt(tsnp), method = 'r', verbose = FALSE))
    df_qc <- .rbind_qc(df_qc, 'TagSNP', tsnp, l_ids)
  }
  l_ids <- lapply(1:2, function(i) match(l_ids[[i]], l_ids_original[[i]]))
  gdata <- genotype_data_subset(gdata, l_ids[[2]], l_ids[[1]])

  list(gdata = gdata, df_qc = df_qc)
}

snprelate_pca <- function(gdata, n_axes = 32, n_cores = 2) {
  stopifnot(inherits(gdata, 'GenotypeData'))

  # get pca
  gds <- request_snpgds_file(gdata)$snpgds
  l_ids <- lapply(paste0('get', c('Scan', 'Snp'), 'ID'), do.call, list(gdata),
    envir = getNamespace('GWASTools'))
  n_axes <- min(c(nscan(gdata), nsnp(gdata)) - 1, n_axes)
  pca <- SNPRelate::snpgdsPCA(gds, l_ids[[1]], l_ids[[2]], num.thread = n_cores,
    eigen.cnt = n_axes, verbose = FALSE)
  if (n_axes == 0) n_axes <- ncol(pca$eigenvect)
  seq_axes <- seq_len(n_axes)
  name_axes <- paste0('PC', seq_axes)

  # get variance
  var <- pca$eigenval[seq_axes] / sum(pca$eigenval, na.rm = TRUE)
  var_percent <- round(var * 100)
  df_variance <- data.frame(matrix(c(var, var_percent), 2, byrow = TRUE),
    PCA_VARNAME = c('Explained_variance', 'Explained_variance_percent'),
    PCA_VARTYPE = 'OTHER')
  names(df_variance)[seq_axes] <- name_axes

  # get obs
  mat_obs <- pca$eigenvect * sqrt(length(l_ids[[1]]))
  df_obs <- data.frame(mat_obs, PCA_VARTYPE = 'OBS',
    PCA_VARNAME = paste0('OBS_', l_ids[[1]]), stringsAsFactors = FALSE)
  names(df_obs)[seq_axes] <- name_axes

  # get vars
  df_vars <- data.frame(PCA_VARNAME = paste0('VAR_', l_ids[[2]]),
    PCA_VARTYPE = 'VAR', stringsAsFactors = FALSE)
  m_vars <- SNPRelate::snpgdsPCASNPLoading(pca, gds, n_cores, FALSE)
  df_vars <- cbind(t(m_vars$snploading), df_vars)
  names(df_vars)[seq_axes] <- name_axes

  df_pca <- df_rbind_all(df_variance, df_obs, df_vars)

  # add obs annotations from gdata
  obs_idxs <- df_pca$PCA_VARTYPE == 'OBS'
  df_annot <- subset(gdata@scanAnnot@data, scanID %in% l_ids[[1]]) 
  for (colname in names(gdata@scanAnnot@data)) {
    df_pca[[colname]] <- NA
    df_pca[obs_idxs, colname] <- as.character(df_annot[[colname]])
    df_pca[[colname]] <- factor(df_pca[[colname]])
  }
  
  class(df_pca) <- c('data.frame', 'pca')

  df_pca
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
  df <- as.data.frame(dplyr::bind_rows(...))
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

###############################################################################
snps_hla <- function(gdata) {
  location <- getSnpVariable(gdata, c('chromosome', 'position'))
  in_hla <- location$chromosome == 6
  in_hla <- in_hla & location$position > 25e6 & location$position < 35e6
}

###############################################################################
.genotype_data_subset <- function(set, gdata) {
  snps_idx <- switch(set,
    hla = which(snps_hla(gdata)),
    non_hla = which(!snps_hla(gdata)),
    getSnpID(gdata))

  if (set == 'no_controls' && !'phenotype' %in% names(gdata@scanAnnot@data)) {
    set <- ''
  }
  scans_idx <- switch(set,
    no_controls = which(gdata@scanAnnot@data$phenotype != 'Control'),
    seq_len(GWASTools::nscan(gdata)))

  genotype_data_subset(gdata, snps_idx, scans_idx)
}
