
.haplo_features <- function(l_haplo, n_cores) {
  snps <- sapply(l_haplo, function(haplo) is.null(dim(haplo)))
  l_haplo[!snps] <- lapply(l_haplo[!snps], function(haplo) {
      haplo_features(haplo)
    })
  l_haplo
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

.merge_haplotypes <- function(m_snps, max_n_snps = 2e2, n_mixtures = 2) {
  # estimate 2 gaussian mixtures, use hclust if too large for performance
  if (ncol(m_snps) > max_n_snps) {
    return(.hclust_classif(hclust(dist(m_snps))))
  }
  clust <- catch_warnings(mclust::Mclust(m_snps, n_mixtures))

  # delete 2 useless warnings
  msgs <- 'best model occurs at the|optimal number of clusters occurs at'
  clust[[2]] <- clust[[2]][!grepl(msgs, clust[[2]])]
  if (length(clust[[2]]) != 0) {
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

  l_peaks <- lapply(l_peaks,
    function(peaks) list(df_snp$probe_id[peaks], peaks))

  l_haplo <- parallel::mclapply(l_peaks, .get_haplos,
    1, 'bedfile', mc.cores = n_cores, mc.set.seed = FALSE)
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
  .system2('shapeit', paste('-B output_plink --seed 1 -T', n_cores, '-O'))
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
