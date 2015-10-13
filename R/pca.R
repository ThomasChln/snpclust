#' Get qb_pca object
#'
#' Create a qb_pca object that stores the pca object and the
#' original data
#'
#' @param data      Data frame subsetted in columns by vars, in rows by obs,
#'                  and passed as parameter data to prcomp.
#' @param id_col    Name of identifier column for observations (must not contain duplicates)
#' @param obs       Row subset of data. No subset if NULL.
#'                    Default: NULL
#' @param vars      Names of columns to be used for the PCA. All the variables
#' which are not specified will be considered as supplementary variables. If NULL
#' all the variables are used except the ID column
#' @param center    Parameter passed to prcomp as center.
#' @param scale     Parameter passed to prcomp as scale..
#'
#' @return List of data, variables, and prcomp object.
#'
#' @export
get_pca    <- function(data, id_col, obs = NULL, vars = NULL, center = TRUE,
  scale = TRUE) {

  check_fx_args(data = '!d+', id_col = '!C1', obs = 'C+', vars = 'C+',
    center = '!B1', scale = '!B1')

  names_data  <- names(data)
  if (is.null(vars)){
    # get all variables if not specified
    vars <- names_data[-match(id_col, names_data)]
  } else {
    # if specified test that vars exist in the dataset
    vars <- unique(vars) # remove potential duplicates
    die_unless(all(c(vars, id_col) %in% names_data),
      "One of the vars or the id_col
        specified is not in the data provided"
    )
  }

  # Test that the vars chosen are all numeric
  nonnum <- which(!sapply(data[, vars], is.numeric))
  die_if(length(nonnum) > 0,
    "Non numeric variable: %s", paste(names(nonnum))
  )

  # Test if unique identifiers
  die_if(any(duplicated(data[, id_col])),
    'Duplicated elements in id_col')


  if (is.null(obs)) {
    # Take all observations
    obs <- data[, id_col]
  } else {
    # take selected observations
    # and cross check if they exist
    obs_missing <- is.na(match(obs, data[, id_col]))
    die_if(any(obs_missing),
      "Could not find the observation: %s", paste(obs[obs_missing])
    )
  }
  #get observation indexes
  obs_idx <- match(obs, data[, id_col])
  data_prcomp <- data[obs_idx, vars]
  rownames(data_prcomp) <- obs
  # Compute pca on active vars and obs
  mypca <- compute_pca_prcomp(data_prcomp,
    center = center,
    scale = scale)

  # create qbpca class
  qb_pca <- new("qb_pca",
    pca = mypca,
    data = data,
    obs = as.character(obs),
    obs_sup = setdiff(as.character(data[,id_col]), obs),
    id_col = id_col,
    vars = vars,
    vars_sup = setdiff(colnames(data), c(id_col, vars))
  )

  qb_pca
}

#' Compute the pca using prcomp
#'
#' @param data data.frame in wide format
#' @param center  Parameter passed to prcomp as center.
#' @param scale   Parameter passed to prcomp as scale..
#' @return prcomp object
#' @author adicara
compute_pca_prcomp <- function(data, center = TRUE, scale = TRUE) {
  prcomp(data,
    retx = TRUE,
    center = center,
    scale. = scale)
}

#' Fortify the original dataset used for the
#' pca with the results of the PCA.
#'
#' @param qb_pca  qb_pca object
#' @param active_obs         Character identifiers of active observations to use. If TRUE no
#'                    subset. If FALSE none are used.
#'                      Default: TRUE
#' @param sup_obs     Character identifiers of suppl. observations to use. If TRUE no
#'                    subset. If FALSE none are used.
#'                      Default: FALSE
#' @param active_vars        Character identifiers of active variables to use. If TRUE no
#'                    subset. If FALSE none are used.
#'                      Default: FALSE
#' @param sup_vars    Character identifiers of suppl. variables to use. If TRUE no
#'                    subset. If FALSE none are used.
#'                      Default: FALSE
#' @param include_pc_variance Logical, includes the explained variance in the data frame
#' @return data.frame of original data with
#' additional columns containing the PCs.
#' The explained variance, variables, and individuals are rbinded in the
#' same data.frame
#' @author adicara
#' @export
pca_fortify <- function(qb_pca,
  active_obs = TRUE,
  sup_obs = TRUE,
  active_vars = TRUE,
  sup_vars = TRUE,
  include_pc_variance = TRUE
) {

  check_fx_args(include_pc_variance = '!B1')
  df_fort <- if (include_pc_variance) {
      get_pca_explained_var(qb_pca)
    } else {
      data.frame()
    }
  df_fort <- df_rbind_all(df_fort,
    get_pca_individuals(qb_pca,
      active = active_obs,
      sup = sup_obs),
    get_pca_vars(qb_pca,
      active = active_vars,
      sup = sup_vars)
  )
  # type the data.frame
  class(df_fort) <- c("qb_pcafort", "data.frame")
  df_fort
}

#' Compute the PC coordinates of supplementary observations
#'
#' @param qb_pca_obj qb_pca object
#' @param obs character vector of obs to retrieve
#' @return data.frame with PC coordinates of sup individuals
#' @author adicara
compute_pca_obs_sup  <- function(qb_pca_obj, obs) {

  #get index of obs_sup. Add rowname to ensure retrieval of ids
  # afterwards
  obs_idx <- match(obs, qb_pca_obj@data[, qb_pca_obj@id_col])
  data_sup <- qb_pca_obj@data[obs_idx, qb_pca_obj@vars]
  rownames(data_sup) <- qb_pca_obj@data[obs_idx, qb_pca_obj@id_col]

  # scale and use rotation matrix then add identifier column
  # using rowname
  if (!is.logical(qb_pca_obj@pca$center)) {
    data_sup <- t(t(data_sup) - qb_pca_obj@pca$center)
  }
  if (!is.logical(qb_pca_obj@pca$scale)) {
    data_sup <- t(t(data_sup) / qb_pca_obj@pca$scale)
  }
  data_sup <- data_sup %*% qb_pca_obj@pca$rotation
  df_out <- as.data.frame(data_sup)
  df_out$PCA_VARNAME <- rownames(df_out)
  df_out$PCA_VARTYPE <- 'OBS_SUP'

  df_out
}

###############################################################################
#' check availability
#'
#' check that elements of selected vector are in the available vector,
#' or if selected is true get all available, if false get none.
#'
#' Jun 5, 2014
#'
#' @param selected logical or character vector, see decription
#' @param available logical or character vector, see decription
#' @param name_error_msg name of the character vector printed in the error message
#' @return selected character vector, see description
#'
#' @author tcharlon
#' @export
check_availability <- function(selected, available, name_error_msg) {
  selected <- if(length(selected) == 1 && selected == TRUE) { # get all
      available
    } else if (length(selected) == 1 && selected == FALSE) { # get none
      character(0)
    } else if (is.character(selected)) { # retrieve selected
      # are they part of the available list ?
      bad_idx <- which(!(selected %in% available))
      die_unless(!length(bad_idx),
        paste("These", name_error_msg, "are not available: %s"),
        paste(selected[bad_idx])
      )
      selected
    } else {
      stop(paste0("Bad specification of ", name_error_msg, ". Use TRUE/FALSE or names"))
    }

  selected
}

#' Extract individuals from the PCA, together
#' with the PCs values
#'
#' @param qb_pca qb_pca object
#' @param active Logical or vector of active observations names to retrieve.
#'              If true all are retrieved, if false none
#' @param sup Logical or vector of supplementary observations names to retrieve.
#'              If true all are retrieved, if false none
#' @return data.frame of individuals
#' @author adicara
#' @export
get_pca_individuals <- function(qb_pca, active = TRUE, sup = TRUE) {

  die_unless(inherits(qb_pca, "qb_pca"),
    "Not a qb_pca object")

  obs <- check_availability(active, qb_pca@obs, 'active observations')
  obs_sup <- check_availability(sup, qb_pca@obs_sup, 'supplementary observations')
  if (!length(obs) && !length(obs_sup)) return(data.frame())

  df_active <- data.frame()
  # get coordinates for active observations
  if (length(obs)) {
#    sqrt_n_obs <- sqrt(length(qb_pca@obs))
    df_active <- as.data.frame(qb_pca@pca$x)# / sqrt_n_obs)
    df_active$PCA_VARNAME <- rownames(df_active)
    df_active$PCA_VARTYPE <- 'OBS'
  }

  df_sup <- data.frame()
  # get coordinates for suppl. observations
  if (length(obs_sup)) {
    df_sup    <- compute_pca_obs_sup(qb_pca, obs_sup)
  }

  df_fortified <- df_rbind_all(df_active, df_sup)

  # merge PC with original data
  df_fortified[[qb_pca@id_col]] <- df_fortified$PCA_VARNAME

  df_fortified <- dplyr::inner_join(qb_pca@data, df_fortified, by = qb_pca@id_col)
  # sort it to get behaviour compatible with dplyr-0.3
  id_col_order <- order(df_fortified[[qb_pca@id_col]])
  df_fortified <- df_fortified[id_col_order, ]

  df_fortified
}

# df_values 2 columns data frame (id, value)
compute_var_coords_numeric <- function(x, df_values) {
  ids <- rownames(x)
  ind <- match(ids, df_values[[1]])
  values <- df_values[ind, 2, drop = TRUE]

  cor(scale(values), x)
}

compute_var_coords_categorical <- function(qb_pca, data_active, varname) {

  values <- data_active[, varname]
  values <- factor(values)
  all_levels <- levels(values)

  has_na <- any(is.na(values))

  # function to compute the coords by categorical variable level
  .compbylevel <- function(var_level) {

    idxget   <- which(values == var_level)
    values <- data_active[idxget, qb_pca@vars]
    new_indiv   <- (colMeans(values) -
        qb_pca@pca$center) / qb_pca@pca$scale
    out <- new_indiv %*% qb_pca@pca$rotation

    #	fix :if contains at least one NA then return NA
    # it could be done before the computation.
    # but in this way we are guaranteed that the
    # cardinality is the same
    if(has_na){
      out[1, ] <- NA
    }
    out <- as.data.frame(out)
    out$PCA_VARNAME <- var_level

    out
  }

  # get coords by levels
  mat_coords  <- plyr::ldply(all_levels, .compbylevel)
  mat_coords$PCA_VARNAME <- paste(varname, mat_coords$PCA_VARNAME, sep = '_|_')
  # integrate variable name in colname
  mat_coords
}

#' Get the coordinates for supplementary variables
#'
#' @param qb_pca qb_pca object
#' @param vars vector of supplementary variable name for which the
#' coordinates need to be computed. By default all the supplementary
#' variables contained in the qb_pca object.
#' @return data.frame containing the PC values for each supplementary variable.
#' Categorical variables are represented by level
#' @author adicara
compute_pca_coords_vars_sup  <- function(qb_pca, vars) {

  check_fx_args(vars = '!C+')

  # Get data for active observations
  data_active <- get_pca_individuals(qb_pca, TRUE, FALSE)

  .compute_sup_var_coords <- function(varname) {
    # Chose computation methofd according to data type
    # numerical or categorical
    df_values <- data_active[, c(qb_pca@id_col, varname), drop = FALSE]

    if (is.numeric(df_values[[2]])) { # numerical
      current_values <- compute_var_coords_numeric(qb_pca@pca$x, df_values)
      current_values <- as.data.frame(current_values)
      current_values$PCA_VARNAME <- varname
    } else {
      current_values <- compute_var_coords_categorical(qb_pca,
        data_active, varname)
    }
    current_values
  }

  # Iterate through supplementary vars and compute coordinates
  dfr_sup_vars <- plyr::ldply(vars, .compute_sup_var_coords)

  dfr_sup_vars
}

#' Extract vars from the PCA, together
#' with the PCs values
#'
#' @param qb_pca qb_pca object
#' @param active Logical or vector of active variable names to retrieve.
#'              If true all are retrieved, if false none
#' @param sup Logical or vector of supplementary variable names to retrieve.
#'              If true all are retrieved, if false none
#' @return data.frame of variables
#' @author adicara
#' @export
get_pca_vars <- function(qb_pca, active = TRUE, sup = TRUE) {

  die_unless(inherits(qb_pca, "qb_pca"),
    "Not a qb_pca object")

  vars <- check_availability(active, qb_pca@vars, 'active variables')
  vars_sup <- check_availability(sup, qb_pca@vars_sup, 'supplementary variables')

  if (!length(vars_sup) & !length(vars)) return(data.frame()) # Nothing to do!

  df_vars_active <- data.frame()
  if (length(vars)) { # get coordinates for active variables
    df_vars_active <- as.data.frame(qb_pca@pca$rotation[vars, ])
    df_vars_active$PCA_VARNAME  <- vars
    df_vars_active$PCA_VARTYPE  <- 'VAR'
  }

  df_vars_sup <- data.frame()
  if(length(vars_sup)) { # get coordinates for active variables
    df_vars_sup <- compute_pca_coords_vars_sup(qb_pca, vars_sup)
    df_vars_sup$PCA_VARTYPE  <- 'VAR_SUP'
  }


  df_out <- df_rbind_all(df_vars_active, df_vars_sup)

  df_out

}


#' Get explained variance for each PC
#'
#' @param qb_pca qb_pca object
#' @return vector of explained variance for each PC expressed in percentage
#' @author adicara
#' @export
get_pca_explained_var <- function (qb_pca) {

  die_unless(inherits(qb_pca, "qb_pca"),
    "Not a qb_pca object")

  explained_var <- qb_pca@pca$sdev ^ 2
  explained_var_percent <- round(100 * explained_var / sum(qb_pca@pca$sdev ^ 2))

  df_explained_var <- as.data.frame(rbind(explained_var, explained_var_percent))
  names(df_explained_var) <- colnames(qb_pca@pca$x)
  df_explained_var$PCA_VARNAME <- c('Explained_variance', 'Explained_variance_percent')
  df_explained_var$PCA_VARTYPE <- 'OTHER'

  df_explained_var
}


# JW fix: Duplicated function in pca_vtest.R
#
## Fit lm between PC and variable of interest
##
## @param dfr data.frame containing variables
## @param pcvecname name of PC column of interest
## @param varvecname name of variable column of interest
## @return data.frame containing t-value and p-value for the
## variable term
## @author adicara
#compute_pca_vtest <- function(dfr, pcvecname, varvecname) {
#
#  check_fx_args(dfr = "!D+", pcvecname = "!C1", varvecname = "!C1")
#
#  frm <- as.formula(paste(pcvecname, "~" ,varvecname))
#  fit <- lm(frm, dfr)
#
#  # get stats of interest
#  dflt <- get_summary(sw_reg(fit, aov=NULL, conf_int = NULL))
#  dflt <- dflt %>%
#    filter((TEST == "t value" | PARAMTYPE == "P_value") &
#        PARAMNAME != "(Intercept)") %>%
#    select(PARAMNAME, PARAMTYPE, TEST, VALUE)
#
#
#  # format the name of the parameter (for cat variables)
#  dflt <- within(dflt,{
#      PARAMTYPE[TEST == "t value"]  <- "T_value"
#      TEST <- NULL
#      .reg <- sprintf("^(%s)(.+)", varvecname)
#      PARAMNAME <- gsub(.reg, "\\1_|_\\2", PARAMNAME)
#    })
#  dflt
#
#}


#' Add supplementary vars to an existing qb_pca object.
#'
#' The new variables need to be supplied as a data.frame
#' containing an identifier column with the same name as
#' the one contained in the qb_pca object.
#' The data is added using a left join (old-data, new-data).
#' This ensures that the original data is not subsetted.
#'
#' @param qb_pca qb_pca object
#' @param df_vars variables data frame
#' @return qb_pca object with additional supplementary variables
#' @author adicara
#' @export
add_pca_sup_vars <- function(qb_pca, df_vars) {

  die_unless(inherits(qb_pca, "qb_pca"),
    "Not a qb_pca object")
  check_fx_args(df_vars = "!D+")

  # Check that identifier column is present in new data
  die_unless(
    df_columns_exist(df_vars, qb_pca@id_col, silent = TRUE),
    "The additonal data must be a data.frame containing the column %s",
    qb_pca@id_col
  )

  # Check that the ids are unique
  die_if(any(duplicated(df_vars[, qb_pca@id_col])),
    "The identifier column %s contains duplicate values",
    qb_pca@id_col
  )

  # list new variables, omitting the ID
  newvars <- colnames(df_vars)[-match(qb_pca@id_col, colnames(df_vars))]

  # extract the data from the pca
  df_original <- qb_pca@data

  # check for duplicated columns
  idxdp <- newvars %in% colnames(df_original)
  die_if(any(idxdp), "Duplicated columns in new data: %s",
    newvars[idxdp])

  # perform left join with new data. This preserves the original data
  df_added <- dplyr::left_join(df_original, df_vars, by = qb_pca@id_col)
  # update qb_pca object
  qb_pca@data <- df_added
  qb_pca@vars_sup <- c(qb_pca@vars_sup, newvars)

  qb_pca

}



setOldClass("prcomp")

#' s4 class to handle PCA parameters and results
#'  
#' @name qb_pca-class
#' @rdname qb_pca-class
#' @exportClass  qb_pca
setClass("qb_pca",
		representation(
				data = "data.frame",
				pca = "prcomp",
				id_col = "character",
				obs = "character",
				obs_sup = "character",
				vars = "character",
				vars_sup = "character"
		)
)
