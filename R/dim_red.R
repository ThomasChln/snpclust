
.haplo_features <- function(l_haplo, n_cores) {
  l_feats <- mclapply(l_haplo, function(haplo) {
      (if (is.null(dim(haplo))) identity else haplo_features)(haplo)
    }, mc.cores = n_cores)
  attributes(l_feats) <- attributes(l_haplo)
  l_feats
}

haplo_features <- function(m_data, order_idxs = FALSE) {
  col_idxs <- lapply(c('hap', 'frequency'), grep, colnames(m_data))
  m_scans <- m_data[, col_idxs[[1]]]
  m_snps <- m_data[, -unlist(col_idxs)]
  merged_haplo <- .merge_haplotypes(m_snps)
  haplos <- as.numeric(factor(.map_haplo_indiv(m_scans, merged_haplo))) - 1
  if (order_idxs) haplos <- list(order(merged_haplo), haplos)

  haplos
}

.map_haplo_indiv <- function(m_scans, haplo) {
  indiv_haplos <- haplo[apply(m_scans, 2, which.max)]
  m_indiv <- matrix(indiv_haplos, ncol = 2)
  m_indiv <- t(apply(m_indiv, 1, sort))
  indiv_haplos <- Reduce(paste0, as.data.frame(m_indiv))
  hap_obs_names <- colnames(m_scans)[grep('hap1code', colnames(m_scans))]
  obs_match <- order(as.numeric(gsub('hap1code.', '', hap_obs_names)))

  indiv_haplos[obs_match]
}

.merge_haplotypes <- function(m_snps, max_n_snps = 2e2, n_mixtures = 2:3) {
  # estimate 2 gaussian mixtures, subset SNPs for performance
  if (ncol(m_snps) > max_n_snps) {
    m_snps <- m_snps[, seq(1, ncol(m_snps), length = max_n_snps)]
  }
  clust <- catch_warnings(mclust::Mclust(m_snps, n_mixtures))
  # delete 2 useless warnings
  msgs <- 'best model occurs at the|optimal number of clusters occurs at'
  clust[[2]] <- clust[[2]][!grepl(msgs, clust[[2]])]
  if (length(clust[[2]])) {
    # bug where mclust gives only one group, use hclust 
    if (length(clust[[2]]) == 1 && grepl('no assignment to', clust[[2]][[1]])) {
      .hclust_classif(hclust(dist(m_snps)))
    } else stop(clust[[2]])
  } else {
    clust[[1]]$classification
  }
}

.hclust_classif <- function(hcl) {
  merged_haplos <- cutree(hcl, 2)
  reorder <- cbind(seq_along(merged_haplos), merged_haplos)[hcl$order, ]
  reorder[, 2] <- as.numeric(factor(reorder[, 2], unique(reorder[, 2])))
  reorder[order(reorder[, 1]), 2]
}

###############################################################################
.estimate_haplotypes <- function(l_peaks, n_cores, df_snp) {
  small_peaks <- sapply(l_peaks, length) < n_cores
  l_peaks <- lapply(l_peaks,
    function(peaks) list(df_snp$probe_id[peaks], peaks))
  l_haplo <- list()
  l_haplo[small_peaks] <- mclapply(l_peaks[small_peaks], .get_haplos,
    1, 'bedfile', mc.cores = n_cores)
  l_haplo[!small_peaks] <- lapply(l_peaks[!small_peaks], .get_haplos,
    n_cores, 'bedfile')
  names(l_haplo) <- names(l_peaks)

  l_haplo
}

.get_haplos <- function(peaks, n_cores, path) {
  ids <- peaks[[1]]
  peaks <- peaks[[2]]
  if (length(peaks) != 1) {
    haplo_mcmc(path, ids, n_cores = n_cores)
  } else {
    peaks
  }
}

###############################################################################
#' haplo_em
#'
#' Get haplotypes with iterative EM from haplo.stats (quadratic complexity).
#' Jan 5, 2015
#'
#' @param gdata       GenotypeData object
#' @param snp_idxs    SNPs indexes
#' @param impute      Function for imputation. If NULL no imputation.
#' @param em_control  haplo.stats::haplo.em.control return
#' @return Numeric matrix (haplotypes x frequency,rs_ids,hap{1,2}code.scan_id)
#' and an optional pvalues column (depending on pvals parameter).
#'
#' @author tcharlon
#' @export
haplo_em <- function(gdata, snp_idxs, impute = sample_impute,
  em_control = haplo.stats::haplo.em.control(min.posterior = 1e-3, n.try = 1)
) {
  stopifnot(inherits(gdata, 'GenotypeData'), inherits(impute, 'function'))

  # Fetch and impute
  gdata <- genotype_data_subset(gdata, snp_idxs)
  m_geno <- t(fetch_genotypes(gdata))
  if (!is.null(impute)) m_geno <- impute(m_geno)

  # Convert to (scans x locus1.allele1,locus1.allele2,...)
  m_geno <- apply(m_geno, 2, function(snp) sapply(c('>', '>='), Reduce, 1, snp))
  m_geno <- matrix(as.numeric(m_geno), GWASTools::nscan(gdata))

  # Compute frequencies and SNP loci
  haplo_em <- haplo.stats::haplo.em(m_geno, GWASTools::getSnpID(gdata), NA,
    control = em_control)
  ord_freqs <- order(haplo_em$hap.prob, decreasing = TRUE)
  frequency <- haplo_em$hap.prob[ord_freqs]
  m_snps <- data.matrix(haplo_em$haplotype[ord_freqs, ])
  m_snps <- abs(m_snps - 1)
  rs_ids <- GWASTools::getSnpVariable(gdata, 'probe_id',
    GWASTools::getSnpID(gdata))
  colnames(m_snps) <- rs_ids

  # Get scans' haplotypes probabilities
  hap_names <- c('hap1code', 'hap2code')
  df_scans <- as.data.frame(haplo_em[c('indx.subj', 'post')])
  l_scans <- lapply(haplo_em[hap_names], .get_em_scans, df_scans)

  m_scans <- t(df_rbind_all(l_scans))

  colnames(m_scans) <- unlist(lapply(hap_names, paste0, '.',
      GWASTools::getScanID(gdata)))
  m_scans <- m_scans[match(ord_freqs, rownames(m_scans)), ]
  m_scans[is.na(m_scans)] <- 0

  # Bind as matrix
  m_haplo <- cbind(frequency, m_snps, m_scans)
  rownames(m_haplo) <- NULL

  m_haplo
}

.get_em_scans <- function(haps, df_scans) {
  df <- cbind(df_scans, haps)
  reshape2::dcast(df, indx.subj ~ haps, value.var = 'post', fill = 0)
}

###############################################################################
#' haplo_mcmc
#'
#' Get haplotypes with MCMC (linear complexity, multithreaded)
#' Jan 5, 2015
#'
#' @param path              BED (optionally tar gz) file path
#' @param probe_ids         Probe ids of SNPs to model
#' @param n_cores           Number of cores to use
#' @return Numeric matrix (haplotypes x frequency,snps,scans)
#'
#' @author tcharlon
#' @export
haplo_mcmc <- function(path, probe_ids, n_cores = 1) {
  if (system('plink', ignore.stdout = TRUE) == 127) stop('plink not found')
  if (system('shapeit') == 127) stop('shapeit not found')

  if (!grepl('^/', path)) path <- file.path(getwd(), path)
  setup_temp_dir()
  # Optional untar
  if (grepl('[.]tar[.]gz$', path)) path <- .untar_bed_targz(path)
  # Subset file with plink
  .system2('plink', paste('--bfile', path, '--snps',
      paste(probe_ids, collapse = ','), '--make-bed --out'))
  # Call shapeit
  .system2('shapeit', paste('-B output_plink -T', n_cores, '-O'))
  # Read output
  df_haplo <- as.data.frame(data.table::fread('output_shapeit.haps'))
  rownames(df_haplo) <- probe_ids
  df_haplo <- df_haplo[-(1:5)]

  .format_mcmc(df_haplo)
}

.format_mcmc <- function(df_haplo) {
  # Compute haplotype frequencies
  n_scan_haplo <- ncol(df_haplo)
  scan_haplos <- Reduce(paste0, as.data.frame(t(df_haplo)))
  scan_haplo_freqs <- table(scan_haplos) / n_scan_haplo

  # Format
  ord_freqs <- order(scan_haplo_freqs, decreasing = TRUE)
  frequency <- scan_haplo_freqs[ord_freqs]
  names_frequency <- names(frequency)
  frequency <- as.vector(frequency)

  # Get haplotypes' SNPs' alleles
  snps <- match(names_frequency, scan_haplos)
  m_snps <- t(df_haplo[, snps])

  # Get scans' most likely haplotypes
  scan_haplos <- match(scan_haplos, names_frequency)
  haplo_id <- paste0('hap', 1:2, 'code')
  scan_id <- sapply(haplo_id, paste0, '.', seq_len(n_scan_haplo / 2))
  scan_id <- as.vector(t(scan_id))
  df_scans <- na.omit(data.frame(scan_haplos, scan_id, 1))
  m_scans <- reshape2::acast(df_scans, scan_haplos ~ scan_id, value.var = 'X1', fill = 0)

  m_haplo <- cbind(frequency, m_snps, m_scans)
  rownames(m_haplo) <- NULL

  m_haplo
}
