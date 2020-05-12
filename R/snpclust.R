
#' snpclust
#'
#' SNPClust performs PCA-based SNP selection and haplotype summarization of
#' nearby correlated SNPS in linear complexity using the SHAPEIT software.
#' Features are then weighted by pricipal component rank and contributions.
#'
#' @param gds       Path of Genomic Data Structure file 
#' @param bed_paths Paths of PLINK binary files.
#'                  If missing, files are generated with the gds file.
#' @param n_axes    Number of principal components to consider
#' @param n_cores   Number of cores to use
#' @param only_pca  Logical, compute only PCA and exit
#' @param keep_tmp  Logical, keep temporary objects
#' @param ...       Passed to snprelate_qc 
#' @return List of slots:
#'           gdata: GenotypeData object of the quality controlled gds,
#'           qc: details about the quality control applied,
#'           pca: PCA applied to the quality controlled gds,
#'           peaks: list of SNP IDs selected in each principal component,
#'           features: matrix of selected SNPs and estimated haplotypes,
#'           features_pca: PCA applied to the features.
#' @export
snpclust <- function(gds, bed_paths, n_axes = 1e2, n_cores = 1,
  only_pca = FALSE, keep_tmp = FALSE, ...) {

  if (missing(bed_paths)) {
    paths <- gds_to_bed(normalizePath(gds))
    on.exit(file.remove(paths))
  }

  gdata <- load_gds_as_genotype_data(gds)
  on.exit(close(gdata))
  if (!keep_tmp) setup_temp_dir()

  snpclust_obj = snprelate_qc(gdata, ...)
  snpclust_obj$pca = snprelate_pca(snpclust_obj$gdata, n_axes, n_cores)

  if (only_pca) return(snpclust_obj)

  haplos <- snpclust_features(snpclust_obj, gdata, paths, n_axes, n_cores)
  snpclust_obj$peaks <- attr(haplos, 'peaks')
  snpclust_obj$max_contributor <- attr(haplos, 'max_contributor')

  snpclust_obj$features <- .haplo_features(haplos, n_cores) 
  snpclust_obj$features <- .haplo_weights(snpclust_obj)

  snpclust_obj$features_qc <- features_qc(snpclust_obj)
  snpclust_obj$features_pca <- get_features_pca(snpclust_obj)

  snpclust_obj
}

# perform quick QC of snpclust features
features_qc <- function(snpclust_obj, weighted = FALSE) {
  m_feats <- as.matrix(snpclust_obj$features)
  polymorphs <- get_polymorphic_cols(m_feats)
  m_feats <- m_feats[, polymorphs]

  m_feats <- if (weighted) {
    t(t(m_feats) * attr(m_feats, 'weights'))
  } else {
    as.matrix(m_feats)
  }

  m_feats <- transitive_tagsnp(sample_impute(m_feats))
}

get_features_pca <- function(snpclust_obj) {
  m_feats <- snpclust_obj$features_qc
  gdata <- snpclust_obj$gdata

  # add observation annotations
  df_scans_annots = gdata_scans_annots(gdata)
  obs_ids <- match(getScanID(gdata), df_scans_annots$scanID)
  df_feats <- cbind.data.frame(m_feats, id = as.character(getScanID(gdata)),
    df_scans_annots[obs_ids, ], stringsAsFactors = FALSE)

  pca_fortify(get_pca(df_feats, 'id', vars = colnames(m_feats)))
}

transitive_tagsnp <- function(m_data, r2 = .8) {
  vars_info <- strsplit(tolower(colnames(m_data)), 'chr_|.loc_|_mb')
  chrs <- sapply(vars_info, '[', 2)
  l_data <- lapply(unique(chrs), function(chr) {
      chr_idx <- which(chrs == chr)
      m_data[, chr_idx, drop = FALSE]
    })
  l_data <- lapply(l_data, .transitive_tagsnp, r2)
  cnames <- unlist(lapply(l_data, colnames))

  matrix(unlist(l_data), ncol = length(cnames),
    dimnames = list(rownames(m_data), cnames))
}

.tagsnp <- function(loc, df_ld, r2) {
  if (all(df_ld[df_ld[[1]] == loc, 3] < r2)) loc
}

.transitive_tagsnp <- function(m_data, r2) {
  col_idx <- if (ncol(m_data) > 1) {
      m_data <- m_data[, ncol(m_data):1]
      m_ld <- stats::cor(m_data) ^ 2
      m_ld[upper.tri(m_ld, TRUE)] <- NA
      df_ld <- stats::na.omit(reshape2::melt(m_ld))
      unlist(lapply(colnames(m_data), .tagsnp, df_ld, r2))
    } else 1

  m_data[, col_idx, drop = FALSE]
}

snpclust_features <- function(snpclust_obj, gdata, paths, n_pcs, n_cores) {
  df_pca = snpclust_obj$pca
  df_vars <- subset(df_pca, DIMRED_VARTYPE == 'VAR')
  snp_ids <- gsub('VAR_', '', df_vars$DIMRED_VARNAME)
  df_snp = gdata_snps_annots(gdata, snp_ids)
  df_snp$probe_id <- as.character(df_snp$probe_id)

  # get peaks
  df_abs_vars <- abs(df_vars[paste0('PC', seq_len(n_pcs))])
  l_peaks <- peak_selection(df_abs_vars, df_snp$chromosome, n_cores)
  l_peaks_copy <- l_peaks
  l_peaks <- unlist(lapply(l_peaks, .separate_peaks, df_snp), FALSE)

  # impute bed, then estimate haplotypes and get SNPs
  df_obs <- subset(df_pca, DIMRED_VARTYPE == 'OBS')
  scan_ids <- as.numeric(gsub('OBS_', '', df_obs$DIMRED_VARNAME))
  tmpfile = tempfile()
  impute_bed(paths, df_snp, l_peaks,
     match(scan_ids, GWASTools::getScanID(gdata)), tmpfile)
  l_haplo <- .estimate_haplotypes(l_peaks, n_cores, df_snp, tmpfile)
  l_haplo <- .add_SNPs(l_haplo, df_vars, gdata, scan_ids)

  # get SNPs in haplotypes and maximum contributions
  pcs <- sapply(strsplit(names(l_peaks), '[.]'), '[', 1)
  contribs <- sapply(seq_along(pcs),
    function(id) df_abs_vars[l_peaks[[id]], pcs[id]])
  attr(l_haplo, 'contribs') <- sapply(contribs, max)
  attr(l_haplo, 'max_contributor') <- sapply(contribs, which.max)
  attr(l_haplo, 'peaks') <- l_peaks_copy
  attr(l_haplo, 'ids') <- lapply(l_peaks, function(peaks) snp_ids[peaks]) 

  l_haplo
}

# contribution and variance weighting
.haplo_weights <- function(snpclust_obj) {
  feats <- snpclust_obj$features
  pca <- snpclust_obj$pca
  pcs <- sapply(strsplit(names(feats), '[.]'), '[', 1)
  pcs_n <- as.numeric(sapply(strsplit(pcs, 'PC'), '[', 2))
  contribs <- attr(feats, 'contribs')
  obs_names <- gsub('OBS_', '', subset(pca, DIMRED_VARTYPE == 'OBS')$DIMRED_VARNAME)
  m_data <- matrix(unlist(feats), length(obs_names),
    dimnames = list(obs_names, names(feats)))
  attributes(m_data) <- c(attributes(m_data), attributes(feats)) 
  names(m_data) <- NULL
  variances <- rep(t(pca[1, unique(pcs_n)]), table(pcs_n))
  attr(m_data, 'weights') <- as.numeric(factor(variances)) * contribs

  m_data
}


.add_SNPs <- function(l_haplo, df_vars, gdata, scan_ids) {
  snps_idxs <- sapply(l_haplo, length) == 1
  snps_ids <- gsub('VAR_', '', df_vars$DIMRED_VARNAME[unlist(l_haplo[snps_idxs])])
  m_geno <- t(fetch_genotypes(gdata))
  m_geno <- m_geno[, match(snps_ids, getSnpID(gdata))]
  m_geno <- m_geno[match(scan_ids, getScanID(gdata)), ]
  l_haplo[snps_idxs] <- as.data.frame(m_geno)
  names(l_haplo)[snps_idxs] <- paste0(names(l_haplo)[snps_idxs], '.', snps_ids)

  l_haplo
}

get_polymorphic_cols <- function(data) {
  cols <- apply(as.matrix(data), 2,
    function(vect) stats::var(stats::na.omit(vect)))
  as.logical(cols)
}
