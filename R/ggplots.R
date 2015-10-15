
###############################################################################
#' ggplot_selection
#'
#' @param l_peaks List returned by peak_selection
#' @param df_pca     QB_pcafort data frame, requires variable contributions
#                    with columns PCA_VARNAME and axes
#' @param axes       Axes to plot
#' @param ncol    ggplot2::facet_wrap
#' @param scales  ggplot2::facet_wrap
#' @param ...     Passed to ggplot2::facet_wrap
#' @param n_points   Maximum number of largest values to plot, for each
#'                   principal component
#' @return ggplot
#' @export
ggplot_selection <- function(l_peaks, df_pca, axes = paste0('PC', 1:10),
  ncol = 5, scales = 'free_y', ..., n_points = 2e3) {

  PCA_VARTYPE <- NULL
  df_vars <- abs(subset(df_pca, PCA_VARTYPE == 'VAR')[axes])
  l_peaks <- l_peaks[names(l_peaks) %in% axes]
  l_peaks_vars <- lapply(seq_along(l_peaks), function(idx) {
      axe <- names(l_peaks)[idx]
      list(axe, l_peaks[[idx]], df_vars[axe])
    })
  l_contribs <- lapply(l_peaks_vars, .ggplot_selection, n_points)
  df_vars <- data.frame(do.call(rbind, l_contribs))
  peak_id_lvls <- names(l_peaks)
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
#' @param gdata      GenotypeData object with snpIDs matching df_vars column
#'                   PCA_VARNAME
#' @param byposition Plot the SNPs by chromosome and poistion or by index
#' @inheritParams ggplot_selection
#' @return ggplot
#' @export
ggplot_manhat <- function(df_pca, gdata, axes = paste0('PC', 1:10),
  byposition = TRUE, n_points = 2e3, ncol = 5, ...) {

  PCA_VARTYPE <- NULL
  df_vars <- subset(df_pca, PCA_VARTYPE == 'VAR')
  ids <- as.numeric(gsub('VAR_', '', df_vars$PCA_VARNAME))
  snpvars <- c('chromosome', 'position', 'probe_id')
  df_snps <- gdata@snpAnnot@data
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
#' @param qb_pcafort  qb_pcafort or qb_pca object (which is fortified)
#' @param groups      Column index or name of qb_pca$data for the grouping of
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
#' @param white_panel Logical, changes standard panel theme
#' @param nas_first   Logical, rearranges data frame to put NAs in groups first
#' @param ... Passed to draw_pca_vars
#' @return            ggplot object
#'
#' @author tcharlon
#' @export
ggplot_pca    <- function(qb_pcafort,
  groups        = NULL,
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
  white_panel   = TRUE,
  nas_first     = TRUE,
  ...) {

  check_fx_args(
    groups = 'C1',
    axes = '!N2',
    ellipses_ci = '!N1',
    ellipses = '!B1',
    label = "C1",
    alpha_fill = '!B1',
    white_panel = '!B1')

  # init those 2 variables to avoid 'NOTE' in the R CMD CHECK
  PCA_VARNAME <- NULL

  # dispatch on qb_pca and qb_pcafort
  if (is(qb_pcafort, 'qb_pca')) {
    qb_pcafort <- pca_fortify(qb_pcafort, obs, obs_sup, vars, vars_sup,
      pc_variance)
  } else if (!is(qb_pcafort, 'qb_pcafort')) {
    stop('qb_pcafort must be either a qb_pca or a qb_pcafort object')
  }

  # get names of axes
  names_axe <- paste0("PC", axes)

  #Format axe labels and add PC explained variance
  axis_labs <- create_pca_axis_labels(qb_pcafort = qb_pcafort, pc_variance = pc_variance,
    axes = axes)

  # Define main data mapping
  mapping <- aes_string(x = names_axe[1],
    y = names_axe[2])

  # Generate plot base
  plt <- ggplot(mapping = mapping) +
    xlab(axis_labs[1]) +
    ylab(axis_labs[2]) +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    geom_vline(aes(xintercept = 0), linetype = "dashed")

  # Display observations
  # check availability of all obs
all_obs <- c(
  check_availability(obs,
    qb_pcafort[qb_pcafort$PCA_VARTYPE == "OBS", 'PCA_VARNAME'],
    'active observations'),
  check_availability(obs_sup,
    qb_pcafort[qb_pcafort$PCA_VARTYPE == "OBS_SUP", 'PCA_VARNAME'],
    'suppl. observations'))

  if (length(all_obs)) {
    df_all_obs <- dplyr::filter(qb_pcafort, PCA_VARNAME %in% all_obs)
    aes_param <- list(group = groups, color = groups, shape = 'PCA_VARTYPE')
    if (alpha_fill) {
      aes_param$alpha <- 'PCA_VARTYPE'
      plt <- plt + scale_alpha_manual(name = 'Observation',
        labels = c(OBS = 'Active', OBS_SUP = 'Suppl.'), values = c(1, 0.5))
    }
    df_all_obs <- reorder_pcafort(df_all_obs = df_all_obs, groups = groups,
      nas_first = nas_first)
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
    if (is(obs_sup, 'logical') && !obs_sup) {
      plt <- plt + guides(alpha = FALSE, shape = FALSE)
    }

    if (!is.null(groups) && !is.numeric(df_all_obs[, groups])) {
      df_all_obs <- df_all_obs[!is.na(df_all_obs[[groups]]), ]
      if (ellipses) {
        plt <- plt + stat_ellipse(aes_string(fill = groups, color = groups),
          data = df_all_obs, ci = ellipses_ci)
      }
      if (link) {
        plt <- plt + geom_line(data = df_all_obs,
          aes_string(group = groups, color = groups)
        )
      }
    }


  }
  # check availability of all vars
all_vars <- c(
  check_availability(vars,
    qb_pcafort[qb_pcafort$PCA_VARTYPE == "VAR", 'PCA_VARNAME'],
    'active variables'),
  check_availability(vars_sup,
    qb_pcafort[qb_pcafort$PCA_VARTYPE == "VAR_SUP", 'PCA_VARNAME'],
    'suppl. variables'))

  # Add variables if requested
  if (length(all_vars)) {
    df_all_vars <- qb_pcafort[match(c('Explained_variance', all_vars), qb_pcafort$PCA_VARNAME), ]
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


  if (white_panel) {
    plt <- plt + theme(panel.background = element_blank(),
      panel.grid.major = element_line(colour = "grey", size = 0.3),
      panel.grid.minor = element_line(colour = "grey", size = 0.2, linetype = 2))
  }

  plt + theme_bw()
}

###############################################################################
#' Get list of pc plots
#'
#' Jul 9, 2014
#'
#' @param pcafort qb_pca or qb_pcafort object
#' @param axes pcs to combine
#' @param max_vars max. number of vars to display
#' @param ... Passed to ggplot_pca
#' @inheritParams ggplot_pca
#' @return list of pca ggplots
#'
#' @author tcharlon
get_pca_panels <- function(pcafort, axes = 1:4, groups = NULL, max_vars = 0L,
  ...) {

  check_fx_args(axes = '!I+', max_vars = '!I1')

  stopifnot(length(axes) > 2)

  # dispatch on qb_pca and qb_pcafort
  if (inherits(pcafort, 'qb_pca')) {
    pcafort <- pca_fortify(pcafort)
  } else if (!inherits(pcafort, 'qb_pcafort')) {
    stop('qb_pcafort must be either a qb_pca or a qb_pcafort object')
  }

  # Get PC plots
  df_axes <- combn(axes, 2)
  vars_idxs <- which(pcafort$PCA_VARTYPE == 'VAR')
  seq_vars <- seq_len(min(max_vars, length(vars_idxs)))
  plots <- apply(df_axes, 2, function(axes) {
      axes <- as.vector(axes)
      pcs <- paste0('PC', axes)
      df_vars <- abs(pcafort[vars_idxs, pcs])
      ord_vars_idxs <- apply(df_vars, 2, order, decreasing = TRUE)[seq_vars, ]
      signif_vars_idxs <- vars_idxs[ord_vars_idxs]
      vars_idxs <- setdiff(vars_idxs, signif_vars_idxs)
      if (length(vars_idxs)) pcafort <- pcafort[-vars_idxs, ]
      ggplot_pca(pcafort, groups, axes, white_panel = FALSE,
          vars = TRUE, scale_sdev = TRUE, ...) +
        labs(x = pcs[1], y = pcs[2]) +
        theme(axis.ticks = element_blank(), axis.text = element_blank())
    })

  # Set title and legend
  n_obs <- sum(pcafort$PCA_VARTYPE == 'OBS')
  n_vars <- sum(pcafort$PCA_VARTYPE == 'VAR')
  plots[[1]] <- plots[[1]] +
    theme(legend.position = 'top') +
    guides(col = guide_legend(nrow = 2))

  groups <- if (!is.null(groups)) paste0('Groups: ', groups, ', ')
  title_text <- paste0(groups, 'Obs.: ', n_obs, ', Var.: ', n_vars)
  plots[[1]] <- plots[[1]] + labs(title = title_text)

  # remove other legends
  plots[-1] <- lapply(plots[-1], '+',
    theme(plot.title = element_blank(), legend.position = 'none'))

  # Put in upper right triangle
  axes_len <- length(axes) - 1
  position_matrix <- matrix(1:(axes_len ^ 2), axes_len)
  triangle_panel <- t(upper.tri(position_matrix, TRUE))
  up_triangle <- which(triangle_panel)
  plots[up_triangle] <- plots
  low_triangle <- !triangle_panel
  plots[which(low_triangle)] <- vector('list', sum(low_triangle))

  plots[t(position_matrix)]
}

###############################################################################
#' Plot combinations of PCs, with variance contributions
#'
#' Jun 22, 2014
#'
#' @param plots       List of ggplots
#' @param only_panels Logical, returns only the plot panels
#' @param legend      Should the legend be kept ?
#' @return Gtable
#'
#' @author tcharlon
grob_pca_panels <- function(plots, only_panels = TRUE, legend = FALSE) {
  stopifnot(inherits(plots, 'list'))

  theme_param <- sapply(1:3, function(i) element_blank())
  names(theme_param) <- c(paste0('axis.', c('text', 'ticks', 'title')))
  if (!legend) theme_param$legend.position <- 'none'
  theme_grob_pca <- do.call('theme', theme_param)
  main_plot <- ggplotGrob(plots[[1]] + theme_grob_pca)
  main_plot <- gtable::gtable_filter(main_plot, 'background|title|guide-box')

  grobs <- NULL
  nulls <- sapply(plots, is.null)
  grobs[nulls] <- lapply(1:sum(nulls), grid::nullGrob)
  sqr_mat_theme <- theme(plot.title = element_blank(), legend.position = 'none',
    panel.grid = element_blank())
  plots[!nulls] <- lapply(plots[!nulls], '+', sqr_mat_theme)
  grobs[!nulls] <- lapply(plots[!nulls], ggplotGrob)
  not_nulls <- which(!nulls)

  # Build grob and put plots
  grid_length <- sqrt(length(grobs))
  unit_mat <- grid::unit(rep(1, grid_length - 1), 'null')
  main_plot <- gtable::gtable_add_rows(main_plot, unit_mat, 3)
  main_plot <- gtable::gtable_add_cols(main_plot, unit_mat, 4)
  main_plot$layout$r[2:(2 + legend)] <- 3 + grid_length
  for (plot_idx in not_nulls) {
    if (only_panels) {
      grobs[[plot_idx]] <- gtable::gtable_filter(grobs[[plot_idx]], 'panel|lab|background')
    }
    main_plot <- gtable::gtable_add_grob(main_plot, grobs[[plot_idx]],
      3 + legend + (plot_idx - 1) %% grid_length,
      4 + (plot_idx - 1) %/% grid_length)
  }

  main_plot
}

#' Plot PCA pairs
#'
#' @inheritParams get_pca_panels
#' @inheritParams ggplot_pca
#' @param ... Passed to get_pca_panels
#' @return NULL
#' @export
plot_pca_pairs <- function(axes, ..., max_vars = 0L, ellipses = TRUE) {
  ggplts <- get_pca_panels(axes = axes, ..., max_vars = max_vars,
    ellipses = ellipses)
  plot(grob_pca_panels(ggplts))
}

reorder_pcafort <- function(df_all_obs, nas_first, groups) {

  if (nas_first && !is.null(groups) && any(is.na(df_all_obs[[groups]]))) {
    nas <- which(is.na(df_all_obs[[groups]]))
    df_all_obs <- rbind(df_all_obs[nas, ], df_all_obs[-nas, ])
  }

  df_all_obs
}

create_pca_axis_labels <- function(qb_pcafort, pc_variance = TRUE, axes,
  names_axe = paste0("PC", axes)) {

  #fix note during checks
  PCA_VARNAME <- NULL

  exp_var   <- dplyr::filter(qb_pcafort, PCA_VARNAME == 'Explained_variance_percent')
  axis_labs <- if (pc_variance && nrow(exp_var)) {
      paste0(names_axe, ' (', exp_var[axes], '% explained var.)')
    } else {
      names_axe
    }

  axis_labs
}

stat_ellipse <- function(mapping = NULL, data = NULL, geom = "polygon",
  position = "identity", alpha = 0.1, ...) {
  # Karl: make sure it is in search path, because even if it is in Depends, if it is
  # loaded via another package, the Depends are not put in search path
  # and proto is need for the StatEllipse class to work
  library(proto)

  ### Stat is not exported. Let's trick R CMD check
  Stat <- get('Stat', envir = getNamespace('ggplot2'))
  GeomPolygon <- get('GeomPolygon', envir = getNamespace('ggplot2'))
  .super <- NULL # dummy

  StatEllipse <- proto(Stat, {
      required_aes <- c("x", "y")
      default_geom <- function(.) GeomPolygon
      objname <- "ellipse"

      calculate_groups <- function(., data, scales, ...){
        .super$calculate_groups(., data, scales,  ...)
      }

      calculate <- function(., data, scales, ci = 0.95, segments = 51, ...) {
        dfn <- 2
        dfd <- length(data$x) - 1
        ellipse <- if (dfd < 3) {
            rbind(c(NA,NA))
          } else {
            v <- cov.wt(cbind(data$x, data$y))
            shape <- v$cov
            center <- v$center
            radius <- sqrt(dfn * qf(ci, dfn, dfd))
            angles <- (0:segments) * 2 * pi / segments
            unit.circle <- cbind(cos(angles), sin(angles))
            t(center + radius * t(unit.circle %*% chol(shape)))
          }

        ellipse <- as.data.frame(ellipse)
        colnames(ellipse) <- c("x","y")
        ellipse
      }
    })

  StatEllipse$new(mapping = mapping, data = data, geom = geom, position = position,
    alpha = alpha,  ...)
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
  colors_vars <- as.numeric(factor(data_vars$PCA_VARTYPE))

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
        label = 'PCA_VARNAME'), data_vars)
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

  #fix note during checks
  PCA_VARTYPE <- NULL

  vars <- grepl('^VAR', df_all_vars$PCA_VARTYPE)

  if (scale_sdev) {
    variance_idx <- match('Explained_variance', df_all_vars$PCA_VARNAME)
    die_if(is.na(variance_idx),
      'scale_sdev is TRUE but no row with PCA_VARTYPE == Explained_variance')
    variance_scale <- unlist(sqrt(df_all_vars[variance_idx, names_axe]))
    df_all_vars[vars, names_axe] <- t(t(df_all_vars[vars, names_axe]) * variance_scale)
  }

  max_var <- max(abs(df_all_vars[vars, names_axe]))
  scale <- .7 * scale / max_var
  df_all_vars[vars, names_axe] <- df_all_vars[vars, names_axe] * scale

  df_all_vars  <- select_(df_all_vars, .dots = c("PCA_VARNAME", "PCA_VARTYPE", names_axe)) %>%
    filter(PCA_VARTYPE != 'OTHER')

  df_all_vars
}



