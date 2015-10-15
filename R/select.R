
peak_selection <- function(df_vars, chromosomes, n_cores = 2) {
  df_abs_vars <- abs(df_vars)
  l_peaks <- lapply(df_abs_vars, .quantile_chromosome_maxs, .6,
    chromosomes)
  lengths <- sapply(l_peaks, length)
  l_vars <- lapply(seq_along(df_abs_vars), function(peak_id) {
      attr(df_abs_vars[[peak_id]], 'max') <- lengths[peak_id]
      df_abs_vars[[peak_id]]
    })
  names(l_vars) <- names(df_abs_vars)
  l_peaks <- mclapply(l_vars, .peaks_with_gmm, mc.cores = n_cores)

  l_peaks[sapply(l_peaks, length) != 0]
}

.separate_peaks <- function(peaks, df_snp) {
  chrs <- df_snp$chromosome[peaks]
  uchrs <- unique(chrs)
  names(uchrs) <- paste0('Chr_', uchrs)
  l_peaks <- lapply(uchrs, function(uchr) {
      .separate_peaks_by_pos(peaks[chrs == uchr], df_snp$position)
    })
  unlist(l_peaks, FALSE)
}

.separate_peaks_by_pos <- function(peaks, position) {
  position <- position[peaks]
  gaps <- which(diff(position) > 1e6)
  gaps <- data.frame(t(cbind(c(1, gaps + 1), c(gaps, length(peaks)))))
  names(gaps) <- paste0('Loc_', round(position[unlist(gaps[1, ])] / 1e6), '_Mb')

  lapply(gaps, function(gap) peaks[seq(gap[1], gap[2])])
}

.peaks_with_gmm <- function(vars, n_vars = min(length(vars), 3e3)) {
  peaks_idxs <- order(vars, decreasing = TRUE)[seq_len(n_vars)]
  for (n_gmm in 2:6) {
    mcl <- Mclust(vars[peaks_idxs], n_gmm)
    class <- which(mcl$uncertainty == 0)
    if (length(class) < attr(vars, 'max')) break
  }

  sort(peaks_idxs[if (n_gmm != 6) class])
}

.quantile_chromosome_maxs <- function(var, quant, chrs) {
  chromosome_maxs <- sapply(unique(chrs), function(chr_id) {
      var_per_chr <- var[chr_id == chrs]
      max(var_per_chr)
    })

  which(var > quantile(chromosome_maxs, quant))
}

.max_peak <- function(peaks, chroms, length_peak) {
  chroms <- chroms[peaks]
  max_peak_chrom <- names(which.max(table(chroms)))
  snps_idxs <- which(chroms == max_peak_chrom)

  if (length(snps_idxs) > length_peak) peaks[snps_idxs]
}

