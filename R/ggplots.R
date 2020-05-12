###############################################################################
#' ggplot_selection
#'
#' Displays the selection of SNPs for each principal component. Peaks is a list
#' of variable indexes corresponding to the pca object.
#'
#' @param peaks    List returned by peak_selection
#' @param pca      PCA data frame, requires variable contributions
#                  with columns DIMRED_VARNAME and axes
#' @param axes     Axes to plot
#' @param ncol     ggplot2::facet_wrap
#' @param scales   ggplot2::facet_wrap
#' @param ...      Passed to ggplot2::facet_wrap
#' @param n_points Maximum number of largest values to plot, for each
#'                 principal component
#' @return ggplot
#' @export
ggplot_selection <- function(peaks, pca, axes = paste0('PC', 1:10),
  ncol = 5, scales = 'free_y', ..., n_points = 2e3) {

  df_vars <- abs(subset(pca, DIMRED_VARTYPE == 'VAR')[axes])
  peaks <- peaks[names(peaks) %in% axes]
  peaks_vars <- lapply(seq_along(peaks), function(idx) {
      axe <- names(peaks)[idx]
      list(axe, peaks[[idx]], df_vars[axe])
    })
  l_contribs <- lapply(peaks_vars, .ggplot_selection, n_points)
  df_vars <- data.frame(do.call(rbind, l_contribs))
  peak_id_lvls <- names(peaks)
  df_vars$peak_id <- factor(df_vars$peak_id, peak_id_lvls)

  ggplot(df_vars, aes_string(x = 'x', y = 'pc', color = 'sel')) + geom_point() +
    facet_wrap(~ peak_id, ncol = ncol, ...) +
    theme(panel.background = element_blank(), panel.grid = element_blank(),
      strip.background = element_blank(), axis.text = element_blank(),
      axis.ticks = element_blank(), legend.position = 'none') +
    labs(x = paste(format(n_points, big.mark =','),
        'most contributing SNPs by index'),
      y = 'Absolute values of contributions')
}

.ggplot_selection <- function(peaks_vars, n_points) {
  df_vars <- peaks_vars[[3]]
  names(df_vars) <- 'pc'
  df_vars$sel <- NA
  df_vars$sel[peaks_vars[[2]]] <- 'Selected'
  pc <- NULL
  threshold <- sort(df_vars$pc, TRUE)[min(nrow(df_vars), n_points)]
  df_vars <- subset(df_vars, pc > threshold)
  seq_vars <- seq_len(nrow(df_vars))

  cbind(df_vars, peak_id = peaks_vars[[1]], x = seq_vars,
    stringsAsFactors = FALSE)
}

###############################################################################
#' ggplot_manhat
#'
#' Displays Manhattan-like plot for SNPs contributions to Principal Components 
#'
#' @param gdata      GenotypeData object with snpIDs matching df_vars column
#'                   DIMRED_VARNAME
#' @param byposition Plot the SNPs by chromosome and poistion or by index
#' @inheritParams ggplot_selection
#' @return ggplot
#' @export
ggplot_manhat <- function(pca, gdata, axes = paste0('PC', 1:10),
  byposition = TRUE, n_points = 2e3, ncol = 5, ...) {

  df_vars <- subset(pca, DIMRED_VARTYPE == 'VAR')
  ids <- as.numeric(gsub('VAR_', '', df_vars$DIMRED_VARNAME))
  snpvars <- c('chromosome', 'position', 'probe_id')
  df_snps <- gdata_snps_annots(gdata)
  df_snps <- df_snps[match(ids, df_snps$snpID), snpvars]
  if (!byposition) df_snps$position <- seq_along(df_snps$position)
  seq_size <- seq_len(min(n_points, nrow(df_vars)))
  df_snps <- plyr::ldply(axes, function(axe) {
      var_contrib <- abs(df_vars[[axe]])
      # select only first with seq_size
      ord_vars <- order(var_contrib, decreasing = TRUE)[seq_size]
      cbind(df_snps[ord_vars, ], axe = axe, var_contrib = var_contrib[ord_vars])
    })
  uniq_chrs <- sort(unique(df_snps$chromosome))
  chrom_lvls <- as.numeric(factor(df_snps$chromosome))
  df_snps$chromosome_color <- factor(chrom_lvls %% 2 + 1)
  if (length(uniq_chrs) > 6) df_snps$position <- round(df_snps$position / 1e6)
  maxs <- c(0, cumsum(sapply(uniq_chrs,
        function(i) max(df_snps$position[df_snps$chromosome == i]))))
  for (chrom_idx in seq_along(uniq_chrs)) {
    snps <- df_snps$chromosome == uniq_chrs[chrom_idx]
    df_snps$position[snps] <- df_snps$position[snps] + maxs[chrom_idx]
  }
  len_maxs <- length(maxs)
  ticks <- diff(maxs) / 2 + maxs[-len_maxs]
  map <- aes_string(x = 'position', y = 'var_contrib',
    color = 'chromosome_color')

  plt <- ggplot(df_snps, map) +
    geom_point() +
    scale_color_discrete(guide = 'none') +
    scale_x_continuous(breaks = ticks, labels = uniq_chrs) +
    facet_wrap(~ axe, ncol = ncol, ...) +
    theme(panel.background = element_blank(), panel.grid = element_blank(),
      strip.background = element_blank(), axis.text.x = element_blank(),
      axis.ticks.x = element_blank(), legend.position = 'none') +
    labs(size = '-log10(p-value)',
      x = paste(format(n_points, big.mark = ','),
        'most contributing SNPs by',
        if (!byposition) 'index' else 'chromosome and position'),
      y = 'Absolute values of contributions')

  plt
}

#' Ggplot a qb_pca object
#'
#' @param pca  pca or qb_pca object (which is fortified)
#' @param group      Column index or name of qb_pca$data for the grouping of
#'                    observations.
#'                      Default: NULL
#' @param axes        Two principal component axes to plot.
#'                      Default: 1:2
#' @param obs         Character identifiers of active observations to use. If TRUE no
#'                    subset. If FALSE none are used.
#'                      Default: TRUE
#' @param obs_sup     Character identifiers of suppl. observations to use. If TRUE no
#'                    subset. If FALSE none are used.
#'                      Default: FALSE
#' @param vars        Character identifiers of active variables to use. If TRUE no
#'                    subset. If FALSE none are used.
#'                      Default: FALSE
#' @param vars_sup    Character identifiers of suppl. variables to use. If TRUE no
#'                    subset. If FALSE none are used.
#'                      Default: FALSE
#' @param ellipses    Should ellipses be plotted ? Not used if groups is NULL.
#'                      Default: FALSE
#' @param ellipses_ci Confidence interval ratio of the ellipses. Not used if
#'                    ellipses is FALSE.
#'                      Default: .95
#' @param label name of column to be used for labelling observation. Can be
#' a supplementary variable or "PCx". By default (NULL), no labels are shown.
#' @param link if TRUE a line is drawn between the observations of the same
#' group. This is useful to identify linked observations
#' @param pc_variance Logical, displays the variance ratio along the axis names
#' @param alpha_fill Logical, should also use alpha to discriminate active and
#'                    additional observations ?
#'                      Default: FALSE
#' @param nas_first   Logical, rearranges data frame to put NAs in groups first
#' @param ... Passed to draw_pca_vars
#' @return            ggplot object
#'
#' @author tcharlon
#' @export
ggplot_pca    <- function(pca,
  group        = NULL,
  axes          = c(1,2),
  obs           = TRUE,
  obs_sup				= FALSE,
  vars          = FALSE,
  vars_sup      = FALSE,
  pc_variance   = TRUE,
  ellipses      = FALSE,
  label         = NULL,
  link 					= FALSE,
  alpha_fill    = FALSE,
  ellipses_ci   = .95,
  nas_first     = TRUE,
  ...) {

  check_fx_args(
    group = 'C1',
    axes = '!N2',
    ellipses_ci = '!N1',
    ellipses = '!B1',
    label = "C1",
    alpha_fill = '!B1')

  # dispatch on qb_pca and pca
  if (methods::is(pca, 'qb_pca')) {
    pca <- pca_fortify(pca, obs, obs_sup, vars, vars_sup,
      pc_variance)
  } else if (!methods::is(pca, 'pca')) {
    stop('pca must be either a qb_pca or a pca object')
  }

  # get names of axes
  names_axe <- paste0("PC", axes)

  #Format axe labels and add PC explained variance
  axis_labs <- create_pca_axis_labels(pca = pca, pc_variance = pc_variance,
    axes = axes)

  # Define main data mapping
  mapping <- aes_string(x = names_axe[1],
    y = names_axe[2])

  # Generate plot base
  plt <- ggplot(mapping = mapping) +
    xlab(axis_labs[1]) +
    ylab(axis_labs[2])

  # Display observations
  # check availability of all obs
  all_obs <- c(
    check_availability(obs,
      pca[pca$DIMRED_VARTYPE == "OBS", 'DIMRED_VARNAME'],
      'active observations'),
    check_availability(obs_sup,
      pca[pca$DIMRED_VARTYPE == "OBS_SUP", 'DIMRED_VARNAME'],
      'suppl. observations'))

  if (length(all_obs)) {
    df_all_obs <- dplyr::filter(pca, DIMRED_VARNAME %in% all_obs)
    aes_param <- list(group = group, color = group, shape = 'DIMRED_VARTYPE')
    if (alpha_fill) {
      aes_param$alpha <- 'DIMRED_VARTYPE'
      plt <- plt + scale_alpha_manual(name = 'Observation',
        labels = c(OBS = 'Active', OBS_SUP = 'Suppl.'), values = c(1, 0.5))
    }
    df_all_obs %<>% reorder_pca(nas_first, group)
    plt <- plt + geom_point(do.call('aes_string', aes_param), df_all_obs) +
      scale_shape_discrete(name = 'Observation',
        labels = c(OBS = 'Active', OBS_SUP = 'Suppl.'))

    # Add ID labels if requested
    if(!is.null(label)) {
      # test that label column  exists
      die_if(!exists(label, where = df_all_obs),
        "Label column %s not found", label)
      # Add geom text
      plt <- plt +
        geom_text(do.call('aes_string',
            c(aes_param, label = label)), df_all_obs)
    }

    #do not show guide if only active observations
    if (methods::is(obs_sup, 'logical') && !obs_sup) {
      plt <- plt + guides(alpha = FALSE, shape = FALSE)
    }

    if (!is.null(group) && !is.numeric(df_all_obs[, group])) {
      df_all_obs <- df_all_obs[!is.na(df_all_obs[[group]]), ]
      if (ellipses) {
        plt <- plt + do.call(ggplot2::stat_ellipse,
            list(mapping = aes_string(color = group),
              data = df_all_obs))
      }
      if (link) {
        plt <- plt + geom_line(data = df_all_obs,
          aes_string(group = group, color = group)
        )
      }
    }
  }
  # check availability of all vars
  all_vars <- c(
    check_availability(vars,
      pca[pca$DIMRED_VARTYPE == "VAR", 'DIMRED_VARNAME'],
      'active variables'),
    check_availability(vars_sup,
      pca[pca$DIMRED_VARTYPE == "VAR_SUP", 'DIMRED_VARNAME'],
      'suppl. variables'))

  # Add variables if requested
  if (length(all_vars)) {
    df_all_vars <- pca[match(c('Explained_variance', all_vars), pca$DIMRED_VARNAME), ]
    # get scale for active vars
    if(obs || obs_sup) {
       circle_sf <-  max(abs(df_all_obs[, names_axe]))
     } else {
       # If no observations are to be shown do not scale variables
       circle_sf <- 1
     }
    var_layer <- draw_pca_vars(df_all_vars, names_axe, circle_sf, ...)
    plt <- plt + var_layer
  }

  plt + theme_bw()
}

reorder_pca <- function(df_all_obs, nas_first, groups) {

  if (nas_first && !is.null(groups) && any(is.na(df_all_obs[[groups]]))) {
    nas <- which(is.na(df_all_obs[[groups]]))
    df_all_obs <- rbind(df_all_obs[nas, ], df_all_obs[-nas, ])
  }

  df_all_obs
}

create_pca_axis_labels <- function(pca, pc_variance = TRUE, axes,
  names_axe = paste0("PC", axes)) {

  exp_var   <- dplyr::filter(pca, DIMRED_VARNAME == 'Explained_variance')
  axis_labs <- if (pc_variance && nrow(exp_var)) {
      paste0(names_axe, ' (', round(exp_var[axes] * 100, 1), '% explained variance)')
    } else {
      names_axe
    }

  axis_labs
}

draw_pca_vars <- function(data_vars, names_axe, scale = 1, circle = FALSE,
  scale_sdev = FALSE) {

  check_fx_args(
    data_vars = '!d+',
    names_axe = '!C2',
    circle = '!B1',
    scale = '!N1')

  data_vars <- scale_pca_vars(df_all_vars = data_vars, names_axe = names_axe,
    scale = scale, scale_sdev = scale_sdev)
  colors_vars <- as.numeric(factor(data_vars$DIMRED_VARTYPE))

  list_geoms <- list(
    geom_segment(aes_string(x = 0, y = 0,
        xend = names_axe[1],
        yend = names_axe[2]),
      data_vars,
      color = c('blue', 'red')[colors_vars],
      arrow = arrow(length = unit(0.1, "inches")),
      lty = 'dashed'
    ),
    geom_text(aes_string(x = names_axe[1], y = names_axe[2],
        label = 'DIMRED_VARNAME'), data_vars)
  )

  if (circle) list_geoms[[3]] <- draw_corr_circle(scale)

  list_geoms
}

draw_corr_circle <- function(scale) {
  angle           <- seq(-pi, pi, length = 500)
  df_circle       <- data.frame(x = sin(angle),
    y = cos(angle)) * scale

  geom_path(aes_string(x = 'x', y = 'y'), df_circle, lty = 'dotted')
}


scale_pca_vars <- function(df_all_vars, names_axe, scale_sdev = FALSE, scale = 1) {

  vars <- grepl('^VAR', df_all_vars$DIMRED_VARTYPE)

  if (scale_sdev) {
    variance_idx <- match('Explained_variance', df_all_vars$DIMRED_VARNAME)
    die_if(is.na(variance_idx),
      'scale_sdev is TRUE but no row with DIMRED_VARTYPE == Explained_variance')
    variance_scale <- unlist(sqrt(df_all_vars[variance_idx, names_axe]))
    df_all_vars[vars, names_axe] <- t(t(df_all_vars[vars, names_axe]) * variance_scale)
  }

  max_var <- max(abs(df_all_vars[vars, names_axe]))
  scale <- .7 * scale / max_var
  df_all_vars[vars, names_axe] <- df_all_vars[vars, names_axe] * scale

  df_all_vars  <- dplyr::select(df_all_vars,
    c("DIMRED_VARNAME", "DIMRED_VARTYPE", names_axe)) %>%
    dplyr::filter(DIMRED_VARTYPE != 'OTHER')

  df_all_vars
}



