
################################################################################
#' snpclust
#'
#' Snpclust performs unsupervised feature selection and summarization based on Principal Component Analysis, Gaussian Mixture Models, and Markov Chain Monte Carlo.
#'
#' @param tar_paths   Path(s) of tar.gz files containing PLINK binary files.
#'                    If missing, files are generated with the gds file.
#' @param gds         Path of Genomic Data Structure file 
#' @param subsets     Character vector for subsetting the dataset.
#'                    '' - no subset
#'                    'hla' - use only SNPs in HLA 
#'                    'non_hla' - use only SNPs outside HLA
#'                    'no_controls' - exclude samples with value Control in
#'                    phenotype 
#' @param n_axes      Number of principal components to consider
#' @param n_cores     Number of cores to use
#' @param ...         Passed to snprelate_qc 
#' @return List of slots: 
#'           pca - PCA applied to the quality controlled dataset,
#'           qc - details about the quality control applied,
#'           gdata - GenotypeData object of the quality controlled dataset,
#'           peaks - list of SNP IDs selected in each principal component,
#'           features - matrix of SNPClust features before weighting,
#'           and features_pca - PCA applied to the weighted SNPClust features.
#'         If several subsets are provided, each slot will be a list with
#'         slots corresponding to each subsets.
#'
#' @author tcharlon
#' @export
snpclust <- function(tar_paths, gds, subsets = '', n_axes = 1e2,
  n_cores = 2, ...) {

  if (missing(tar_paths)) {
    tar_paths <- gds_to_bedtargz(normalizePath(gds))
    on.exit(file.remove(tar_paths))
  }

  # subset gdata, qc and pca
  gdata <- load_gds_as_genotype_data(gds)
  on.exit(close(gdata))
  pca_objs <- lapply(subsets, .snpclust_qc_pca, gdata, n_axes, n_cores, ...)
  snpclust_obj <- list(pca = lapply(pca_objs, '[[', 'pca'),
    qc = lapply(pca_objs, '[[', 'qc'), gdata = lapply(pca_objs, '[[', 'gdata'))

  # haplotypes estimation, merging, and weighting by PC rank and contributions 
  haplos <- lapply(snpclust_obj$pca, .snpclust_features, gdata, tar_paths, 0,
    n_cores)
  snpclust_obj$peaks <- lapply(haplos, attr, 'peaks')
  snpclust_obj$features <- lapply(haplos, .haplo_features, n_cores)
  snpclust_obj$features <- lapply(seq_along(subsets), .haplo_weights,
    snpclust_obj)
  snpclust_obj$features_pca <- lapply(seq_along(subsets), get_features_pca,
    snpclust_obj)

  # unlist slots if only one subsets 
  if (length(subsets) == 1) snpclust_obj <- lapply(snpclust_obj, '[[', 1)

  snpclust_obj
}

.snpclust_qc_pca <- function(subset, gdata, n_axes, n_cores, ...) {
  gdata <- .genotype_data_subset(subset, gdata)
  qced_gdata <- snprelate_qc(gdata, ...)
  pca <- snprelate_pca(qced_gdata$gdata, n_axes, n_cores)
  list(pca = pca, qc = qced_gdata$df_qc, gdata = qced_gdata$gdata)
} 

get_features_pca <- function(idx, snpclust_obj) {
  m_feats <- snpclust_obj$features[[idx]]
  m_feats <- t(t(scale(m_feats)) * attr(m_feats, 'weights'))
  m_feats <- transitive_tagsnp(sample_impute(m_feats))
  obs_ids <- getScanID(snpclust_obj$gdata[[idx]])
  df_feats <- cbind.data.frame(m_feats, id = paste(1:nrow(m_feats)),
    # add observation annotations
    snpclust_obj$gdata[[idx]]@scanAnnot@data[obs_ids, ],
    stringsAsFactors = FALSE)

  get_pca(df_feats, 'id', vars = colnames(m_feats))
}

###############################################################################
transitive_tagsnp <- function(m_data, r2 = .8) {
  vars_info <- strsplit(colnames(m_data), 'Chr_|.Loc_|_Mb')
  chrs <- sapply(vars_info, '[', 2)
  l_data <- lapply(unique(chrs), function(chr) {
      chr_idx <- which(chrs == chr)
      m_data[, chr_idx, drop = FALSE]
    })
  l_data <- lapply(l_data, .transitive_tagsnp, r2)
  cnames <- unlist(lapply(l_data, colnames))

  matrix(unlist(l_data), ncol = length(cnames), dimnames = list(rownames(m_data), cnames))
}

.tagsnp <- function(loc, df_ld, r2) {
  if (all(df_ld[df_ld[[1]] == loc, 3] < r2)) loc
}

.transitive_tagsnp <- function(m_data, r2) {
  col_idx <- if (ncol(m_data) > 1) {
      m_data <- m_data[, ncol(m_data):1]
      m_ld <- cor(m_data) ^ 2
      m_ld[upper.tri(m_ld, TRUE)] <- NA
      df_ld <- na.omit(reshape2::melt(m_ld))
      unlist(lapply(colnames(m_data), .tagsnp, df_ld, r2))
    } else 1

  m_data[, col_idx, drop = FALSE]
}

.snpclust_features <- function(df_pca, gdata, tar_paths, n_pcs, n_cores) {
  PCA_VARTYPE <- NULL
  df_vars <- subset(df_pca, PCA_VARTYPE == 'VAR')
  snp_ids <- gsub('VAR_', '', df_vars$PCA_VARNAME)
  df_snp <- getSnpVariable(gdata, c('chromosome', 'position', 'probe_id'),
    snp_ids)
  df_snp$probe_id <- as.character(df_snp$probe_id)

  # get peaks
  if (!n_pcs) n_pcs <- length(grep('PC[0-9]', colnames(df_vars)))
  df_abs_vars <- abs(df_vars[paste0('PC', seq_len(n_pcs))])
  l_peaks <- peak_selection(df_abs_vars, df_snp$chromosome, n_cores)
  l_peaks_copy <- l_peaks
  l_peaks <- unlist(lapply(l_peaks, .separate_peaks, df_snp), FALSE)

  # Untar and impute bed, then estimate haplotypes and get SNPs
  df_obs <- subset(df_pca, PCA_VARTYPE == 'OBS')
  scan_ids <- as.numeric(gsub('OBS_', '', df_obs$PCA_VARNAME))
  setup_temp_dir()
  .untar_impute_bed(tar_paths, df_snp, l_peaks, scan_ids)
  l_haplo <- .estimate_haplotypes(l_peaks, n_cores, df_snp)
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

# contribution and variance pondering
.haplo_weights <- function(run, tables) {
  feats <- tables$features[[run]]
  pca <- tables$pca[[run]]
  pcs <- sapply(strsplit(names(feats), '[.]'), '[', 1)
  pcs_n <- as.numeric(sapply(strsplit(pcs, 'PC'), '[', 2))
  contribs <- attr(feats, 'contribs')
  PCA_VARTYPE <- NULL
  obs_names <- gsub('OBS_', '', subset(pca, PCA_VARTYPE == 'OBS')$PCA_VARNAME)
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
  snps_ids <- gsub('VAR_', '', df_vars$PCA_VARNAME[unlist(l_haplo[snps_idxs])])
  m_geno <- t(fetch_genotypes(gdata))
  m_geno <- m_geno[, match(snps_ids, getSnpID(gdata))]
  m_geno <- m_geno[match(scan_ids, getScanID(gdata)), ]
  l_haplo[snps_idxs] <- as.data.frame(m_geno)
  names(l_haplo)[snps_idxs] <- paste0(names(l_haplo)[snps_idxs], '.', snps_ids)

  l_haplo
}
