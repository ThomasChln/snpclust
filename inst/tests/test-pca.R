context("PCA")
suppressPackageStartupMessages(library(magrittr))

# Build test dataset
.pca_test_set <- function () {
  cbind(iris,
      SPECIES_NUMERIC  = as.numeric(as.factor(iris$Species)),
      SAMPLEID = paste0("ID", seq_len(nrow(iris))), stringsAsFactors = FALSE)
}



# Build pca test set
.pca_generate <- function (obs = NULL, vars = NULL) {

  dt <- .pca_test_set()

  if(is.null(vars)) {
    vars <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
  }

  if(is.null(obs)) {
    obs <- as.character(dt$SAMPLEID)
  }

  pca <- get_pca(dt,
      obs = obs,
      vars = vars,
      id_col = "SAMPLEID"
  )


  pca
}

.generate_sup_var_names <- function(df_sups, separator = "_|_") {

  names <- sapply(names(df_sups), function(nm){
        dtr <- df_sups[,nm]
        if(is.numeric(dtr)) {
          return(nm)
        } else {
          dtr <- factor(dtr)
          mynames <- paste(nm, levels(dtr), sep = separator)
          return(mynames)
        }
      })

  unname(unlist(names))
}


# Test the pca procedure and generation of qb_pca object
.test_get_pca  <- function() {

  mydata <- .pca_test_set()

  # Normal run =======
 vars  <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
  mypca <- get_pca(mydata,
    vars = vars,
    id_col = "SAMPLEID"
  )

  # Variables are kept
  expect_identical(mypca@vars, vars)

  # Matches with results of compute pca
  pr_manual <- snpclust:::compute_pca_prcomp(mydata[,vars])
    expect_true(all.equal(
                mypca@pca$x,
                pr_manual$x,
                check.attributes = FALSE))

  ## Specify less variables =====
  vars_less  <- c("Sepal.Length", "Sepal.Width")
  mypcalessvars <- get_pca(mydata,
      vars = vars_less,
      id_col = "SAMPLEID"
  )
  # Variables are kept
  expect_identical(mypcalessvars@vars, vars_less)

  # Matches with results of compute pca
  pr_manual <- snpclust:::compute_pca_prcomp(mydata[,vars_less])
  expect_true(all.equal(
          mypcalessvars@pca$x,
          pr_manual$x,
          check.attributes = FALSE))


  # Bad variables or id col ======
  expect_error(
      get_pca(mydata,
          vars = c(vars, "BANANA"),
          id_col = "SAMPLEID"),
      "One of the vars"
  )

  expect_error(
      get_pca(mydata,
          vars = c(vars),
          id_col = "BANANA"),
      "One of the vars"
  )

  # Non numeric variables =====
  expect_error(
      get_pca(mydata,
          vars = c(vars, "Species"),
          id_col = "SAMPLEID"),
      "Non numeric variable: Species"
  )

  # subset observations =====
  myobs <- paste0("ID", c(1:10, 20:30))
  myactive_dataset <- with(mydata, mydata[SAMPLEID %in% myobs, vars])
  mypcasubsetmanual <- snpclust:::compute_pca_prcomp(myactive_dataset)

  mypca_subset  <- get_pca(mydata,
                          vars = vars,
                          obs = myobs,
                          id_col = "SAMPLEID")

  # obs and obs_sup
  expect_equal(mypca_subset@obs, myobs)
  expect_equal(mypca_subset@obs_sup, with(mydata, as.character(
                                      SAMPLEID[!SAMPLEID %in% myobs])))
  # test subsetted matrix
  expect_true(
      all.equal(
          mypca_subset@pca$x,
          mypcasubsetmanual$x,
          check.attributes = FALSE
          )
      )

    # bad subset observations =====
   myobsbad <- paste0("ID", c(1:10, 20:30, "_BANANA"))

   expect_error(
     get_pca(mydata,
         vars = vars,
         obs = myobsbad,
         id_col = "SAMPLEID"),
          "Could not find the observation"
    )

    # sup vars has a NA
    mydatahasna <- mydata
    mydatahasna[1, "Species"]  <- NA

    pcanasupvar <- get_pca(mydatahasna,
                           vars = vars,
                           id_col = "SAMPLEID")

    expect_true(all.equal(
           pcanasupvar@pca$x,
           mypca@pca$x,
           check.attributes = FALSE))

}

# Test retrieval of indivs
.test_get_pca_individuals <- function() {


  pcares  <- .pca_generate()

  df_fort <- get_pca_individuals(pcares)
  expect_true(
      unique(df_fort$PCA_VARTYPE) == 'OBS'
  )
  expect_equal(
      dplyr::filter(df_fort, PCA_VARTYPE == 'OBS') %>% nrow(),
      150L
  )
  #are all the variables there
  original_data_colnames <- colnames(pcares@data)
  expect_true(
    all(original_data_colnames %in% colnames(df_fort))
  )


  # try with variable subsets
  myobs <- paste0("ID", c(1:10, 20:30))
  pcaressub  <- .pca_generate(myobs)
  df_fort_sub <- get_pca_individuals(pcaressub)
  expect_identical(
      unique(df_fort_sub$PCA_VARTYPE), c('OBS', 'OBS_SUP')
  )
  # test counts obs obs_sup
  expect_true(
      all(table(df_fort_sub$PCA_VARTYPE) == c(length(myobs), 150-length(myobs)))
  )
  expect_equal(
      dplyr::filter(df_fort, PCA_VARTYPE == 'OBS') %>% nrow(),
      150L
  )

  # Ask not to compute supp indivs
  df_fort_sub_nosup <- get_pca_individuals(pcaressub, sup = FALSE)
  expect_true(
      unique(df_fort_sub_nosup$PCA_VARTYPE) == 'OBS'
  )
  expect_equal(
      dplyr::filter(df_fort_sub_nosup, PCA_VARTYPE == 'OBS') %>% nrow(),
      length(myobs)
  )


  # shuffle data
  modified_pca_res <- pcares
  modified_pca_res@data <- modified_pca_res@data[sample(1:150, 150), ]

  fort_orig <- pca_fortify(pcares, active_vars = FALSE, sup_vars = TRUE) %>%
                dplyr::arrange(SAMPLEID)

  fort_modif <- pca_fortify(modified_pca_res, active_vars = FALSE, sup_vars = TRUE, sup_obs = FALSE) %>%
                dplyr::arrange(SAMPLEID)
  expect_identical(fort_orig, fort_modif)

  # shuffle pca data.frame
  modified_pca_res <- pcares


  modified_pca_res@pca$x <- modified_pca_res@pca$x[sample(1:150, 150), ]

  # added by karl, otherwise identical(fort_orig, fort_modif2) failed in dplyr 0.3 and 0.4
  fort_orig2 <- pca_fortify(modified_pca_res, active_vars = FALSE, sup_vars = TRUE) %>%
    dplyr::arrange(SAMPLEID)

  fort_modif2 <- pca_fortify(modified_pca_res, active_vars = FALSE, sup_vars = TRUE, sup_obs = FALSE) %>%
                dplyr::arrange(SAMPLEID)
  expect_identical(fort_orig2, fort_modif2)

  # shuffle obs
  modified_pca_res <- pcares
  modified_pca_res@obs <- modified_pca_res@obs[sample(1:150, 150)]
  fort_modif <- pca_fortify(modified_pca_res, active_vars = FALSE, sup_vars = TRUE, sup_obs = FALSE) %>%
      dplyr::arrange(SAMPLEID)
  expect_identical(fort_orig, fort_modif)

  # Perform shuffle with indiv sups
  modified_pca_res_sub <- pcaressub
  fort_orig_sub <- pca_fortify(modified_pca_res_sub,
                        active_vars = FALSE, sup_vars = TRUE) %>%
                        dplyr::arrange(SAMPLEID)

  nrowssp <- nrow(modified_pca_res_sub@pca$x)
  modified_pca_res_sub@pca$x <- modified_pca_res_sub@pca$x[sample(1:nrowssp, nrowssp), ]
  fort_modif_sub <- pca_fortify(modified_pca_res_sub, active_vars = FALSE,
                              sup_vars = TRUE) %>%
                dplyr::arrange(SAMPLEID)
  expect_identical(fort_orig_sub, fort_modif_sub)


}

# Test retrieval of vars
.test_get_pca_vars <- function() {

  myvars <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
  pcares  <- .pca_generate(vars = myvars)
  df_fort <- get_pca_vars(pcares, active = TRUE, sup = FALSE)

  expect_identical(df_fort$PCA_VARNAME, myvars)
  expect_true(all(df_fort$PCA_VARTYPE == 'VAR'))

  # less active variables
  myvarssub <- c("Sepal.Length", "Sepal.Width", "Petal.Width")
  pcares  <- .pca_generate(vars = myvarssub)
  df_fort <- get_pca_vars(pcares, sup = FALSE)

  expect_identical(df_fort$PCA_VARNAME, myvarssub)
  expect_true(all(df_fort$PCA_VARTYPE == 'VAR'))

  # perform pca on all variables but retrieve only some
  df_fort <- get_pca_vars(pcares, active = myvarssub, sup = FALSE)
  expect_identical(df_fort$PCA_VARNAME, myvarssub)
  expect_true(all(df_fort$PCA_VARTYPE == 'VAR'))

  # retieve also supplementary variables
  df_fort <- get_pca_vars(pcares, active = myvarssub, sup = TRUE)
  # active vars
  df_vort_VAR <- df_fort %>% dplyr::filter(PCA_VARTYPE == "VAR")
  expect_identical(df_vort_VAR$PCA_VARNAME, myvarssub)
  #name of sup vars (! some are categorical)
  df_vort_VARSUP <- df_fort %>% dplyr::filter(PCA_VARTYPE == "VAR_SUP")
  mysupvarnames <- .generate_sup_var_names(pcares@data[, pcares@vars_sup])
  expect_equal(df_vort_VARSUP$PCA_VARNAME, mysupvarnames)

  # retieve only supplementary variables
  df_fort <- get_pca_vars(pcares, active = FALSE, sup = TRUE)
  expect_true(all(df_vort_VARSUP$PCA_VARTYPE == "VAR_SUP"))

  # retieve selected sup variables
  df_fort <- get_pca_vars(pcares, active = FALSE, sup = pcares@vars_sup[1])
  expect_true(all(df_vort_VARSUP$PCA_VARTYPE == "VAR_SUP"))
  expect_true(all(df_vort_VARSUP$PCA_VARTNAME == pcares@vars_sup[1]))


}

.compute_pca_coords_vars_sup <- function() {

  mypca <- .pca_generate()
  df_sup_num <- snpclust:::compute_pca_coords_vars_sup(mypca, "SPECIES_NUMERIC")
  expect_true(nrow(df_sup_num) == 1)
  expect_equal(df_sup_num$PCA_VARNAME, "SPECIES_NUMERIC")

  df_sup_cat <- snpclust:::compute_pca_coords_vars_sup(mypca, "Species")
  expect_true(nrow(df_sup_cat) == 3)
  nms <- .generate_sup_var_names(mypca@data[,"Species", drop = FALSE])
  expect_equal(df_sup_cat$PCA_VARNAME, nms[, 1, drop=TRUE])

  # check the actual calculations for numeric and first two lines of categorical
  df_sup <- rbind(df_sup_num, df_sup_cat[1:2, ])
  expected <- structure(list(PC1 = c(0.944665176870244, -2.21732491513681,
        0.494790440357863), PC2 = c(0.0118216838874303, -0.2879627489894,
        0.548333521629207), PC3 = c(-0.193748093327065, 0.0426960630417206,
        0.0958085424866922), PC4 = c(-0.0179874841448261, 0.0182795184349202,
        -0.0302387861298804), PCA_VARNAME = c("SPECIES_NUMERIC", "Species_|_setosa",
        "Species_|_versicolor")), .Names = c("PC1", "PC2", "PC3", "PC4",
      "PCA_VARNAME"), row.names = c(NA, 3L), class = "data.frame")

  expect_equal(df_sup, expected)


  # sup vars has a NA in categorical (return NA)
  mydatahasna <- .pca_test_set()
  mydatahasna[1, "Species"]  <- NA
  pcanasupvar <- get_pca(mydatahasna,
                        vars = c("Sepal.Length", "Sepal.Width",
                                 "Petal.Length", "Petal.Width"),
                        id_col = "SAMPLEID")
  df_sup_catmiss <- snpclust:::compute_pca_coords_vars_sup(pcanasupvar, "Species")
  expect_true(all(is.na(df_sup_catmiss[, 1:3])))

  # sup vars has a NA in numeric (return NA)
  mydatahasnanum <- .pca_test_set()
  mydatahasnanum[1, "Petal.Width"]  <- NA
  pcanasupvarnum <- get_pca(mydatahasnanum,
      vars = c("Sepal.Length", "Sepal.Width",
                "Petal.Length"),
      id_col = "SAMPLEID")
  df_sup_nummiss <- snpclust:::compute_pca_coords_vars_sup(pcanasupvarnum,
                                                        "Petal.Width")
  expect_true(all(is.na(df_sup_nummiss[, 1:3])))
}

.test_pca_fortify <- function() {

  # test fortify agains manual rbind

  mypca <- .pca_generate()
  df_fort <- pca_fortify(mypca)
  manual_fortify <- plyr::rbind.fill(get_pca_explained_var(mypca),
                                get_pca_individuals(mypca),
                                get_pca_vars(mypca,
                                            active = TRUE, sup = TRUE))
  expect_equivalent(df_fort, manual_fortify)

  mypca <- .pca_generate()
  df_fort <- pca_fortify(mypca, sup_obs = FALSE)
  manual_fortify <- plyr::rbind.fill(get_pca_explained_var(mypca),
    get_pca_individuals(mypca,
          sup = FALSE),
      get_pca_vars(mypca,
          active = TRUE, sup = TRUE))
  expect_equivalent(df_fort, manual_fortify)

  mypca <- .pca_generate()
  df_fort <- pca_fortify(mypca, sup_obs = FALSE, active_vars = FALSE)
  manual_fortify <- plyr::rbind.fill(get_pca_explained_var(mypca),
    get_pca_individuals(mypca,
          sup = FALSE),
      get_pca_vars(mypca,
          active = FALSE, sup = TRUE))
  expect_equivalent(df_fort, manual_fortify)

  mypca <- .pca_generate()
  df_fort <- pca_fortify(mypca, sup_obs = FALSE, active_vars = TRUE,
                          sup_vars = FALSE)
  manual_fortify <- plyr::rbind.fill(get_pca_explained_var(mypca),
    get_pca_individuals(mypca,
          sup = FALSE),
      get_pca_vars(mypca,
          active = TRUE, sup = FALSE))
  expect_equivalent(df_fort, manual_fortify)

  # no variables at all
  df_fort <- pca_fortify(mypca, sup_obs = FALSE, active_vars = FALSE,
    sup_vars = FALSE)
  manual_fortify <- plyr::rbind.fill(get_pca_explained_var(mypca),
    get_pca_individuals(mypca,
      sup = FALSE),
    get_pca_vars(mypca,
      active = FALSE, sup = FALSE))

  # no variance
  df_fort <- pca_fortify(mypca, include_pc_variance = FALSE)
  manual_fortify <- plyr::rbind.fill(get_pca_individuals(mypca),
    get_pca_vars(mypca))

  expect_equivalent(df_fort, manual_fortify)
  #bad object
  expect_error(
    pca_fortify(iris),
    "Not a qb_pca object"
  )

}





.test_get_pca_explained_var <- function() {

  dt <- .pca_test_set()
  vars <- colnames(dt)[1:4]
  pcares  <- get_pca(dt,
    vars = vars,
    id_col = "SAMPLEID"
  )

  expect_equal(
    get_pca_explained_var(pcares),
    structure(list(PC1 = c(2.918497816532, 73), PC2 = c(0.91403047146807,
          23), PC3 = c(0.146756875571315, 4
        ), PC4 = c(0.0207148364286192, 1), PCA_VARNAME = c("Explained_variance",
          "Explained_variance_percent"), PCA_VARTYPE = c("OTHER", "OTHER"
        )), .Names = c("PC1", "PC2", "PC3", "PC4", "PCA_VARNAME", "PCA_VARTYPE"
      ), row.names = c("explained_var", "explained_var_percent"), class = "data.frame")
  )

}

.test_add_pca_sup_vars <- function() {

  # get the pca
  mypca <- .pca_generate()

  # generate variables to add
  dforig_data <- pca_fortify(mypca, include_pc_variance = FALSE,
                                    active_vars = FALSE,
                                    sup_vars = FALSE)

  dfmodified <- within(dforig_data, {
        NEW_SAMPLEID <- SAMPLEID
        NEW_SPECIES <- Species
      })
  select_new_data_cols <- c("SAMPLEID", "NEW_SAMPLEID", "NEW_SPECIES")
  df_newdata <- dfmodified[, select_new_data_cols]

  # update the pca
  myupdatedpca <- add_pca_sup_vars(mypca, df_newdata)

  # test if new fortified pca matches with the one previosly generated
  myupdatedpca_fortify <- pca_fortify(myupdatedpca, include_pc_variance = FALSE,
                            active_vars = FALSE,
                            sup_vars = FALSE)

  expect_identical(myupdatedpca_fortify[colnames(dfmodified)],
                    dfmodified)

  # vars_sup updated
  expect_identical(
      myupdatedpca@vars_sup, c(mypca@vars_sup, c("NEW_SAMPLEID", "NEW_SPECIES"))
    )

  # bad id columns
  df_newdata_bad_idcol <- df_newdata
  colnames(df_newdata_bad_idcol)[1] <- "BANANAID"
  expect_error(add_pca_sup_vars(mypca, df_newdata_bad_idcol),
        "containing the column SAMPLEID"
        )

  # non-unique IDs
  df_newdata_dupids <- df_newdata
  df_newdata_dupids[1:2, "SAMPLEID"] <- "ID1"
  expect_error(add_pca_sup_vars(mypca, df_newdata_dupids),
      "The identifier column SAMPLEID"
  )

  #duplicated columns
  df_newdata_dupcols <- df_newdata
  df_newdata_dupcols$Sepal.Length <- TRUE
  expect_error(add_pca_sup_vars(mypca, df_newdata_dupcols),
      "Duplicated columns in new data: Sepal.Length"
  )

  # Test the left join =========>

  ## some new data missing 10 observations
  df_newdata_miss <- df_newdata[-c(1:10), ]
  res_missing <- add_pca_sup_vars(mypca, df_newdata_miss)
  missing_fortify <- pca_fortify(res_missing, sup_obs = TRUE, include_pc_variance = FALSE,
    active_vars = FALSE,
      sup_vars = FALSE)
  # orinal data rows are mantained
  expect_true(nrow(missing_fortify) == nrow(dfmodified))
  # missing values are NA
  expect_equal(
   length(which(is.na(missing_fortify$NEW_SAMPLEID))), 10
  )
  expect_equal(
      length(which(is.na(missing_fortify$NEW_SPECIES))), 10
  )

  ## some new data with additional observations
  # they should not be added
  addobs <- df_newdata[1:5, ]
  addobs$SAMPLEID <- paste0("XXX", 1:5)
  df_newdata_toomuch <- rbind(df_newdata, addobs)
  res_toomuch <- add_pca_sup_vars(mypca, df_newdata_toomuch)
  toomuch_fortify <- pca_fortify(res_toomuch, sup_obs = TRUE,
    include_pc_variance = FALSE,
      active_vars = FALSE,
      sup_vars = FALSE)
  # orinal data rows are mantained
  expect_true(nrow(toomuch_fortify) == nrow(dfmodified))
  expect_identical(toomuch_fortify$NEW_SAMPLEID,
                   toomuch_fortify$SAMPLEID
                  )

}

.test_compute_pca_obs_sup <- function() {
  # check first 2 lines of result with 20 random supplementary observations
  set.seed(1)
  all_obs <- as.character(.pca_test_set()$SAMPLEID)
  actives <- sample(all_obs, 130)
  sups <- setdiff(all_obs, actives)
  pca_iris_supp_obs <- .pca_generate(
    obs = actives)

  df_sup_obs <- snpclust:::compute_pca_obs_sup(pca_iris_supp_obs, sups)[1:2, ]

  expected  <- structure(list(PC1 = c(-2.1000241837692, -2.13252150474489),
      PC2 = c(-1.1274375877109, -1.56924447981553), PC3 = c(-0.270340238555629,
        -0.0118193253118264), PC4 = c(0.016139736047961, 0.195336227514343
      ), PCA_VARNAME = c("ID11", "ID17"), PCA_VARTYPE = c("OBS_SUP",
        "OBS_SUP")), .Names = c("PC1", "PC2", "PC3", "PC4", "PCA_VARNAME",
      "PCA_VARTYPE"), row.names = c("ID11", "ID17"), class = "data.frame")
  expect_equal(df_sup_obs, expected)
}

test_that(".test_get_pca", .test_get_pca())
test_that(".test_get_pca_individuals", .test_get_pca_individuals())
test_that(".test_get_pca_vars", .test_get_pca_vars())
test_that(".test_get_pca_explained_var", .test_get_pca_explained_var())
test_that(".compute_pca_coords_vars_sup", .compute_pca_coords_vars_sup())
test_that(".test_pca_fortify", .test_pca_fortify())
test_that(".test_add_pca_sup_vars", .test_add_pca_sup_vars())
test_that(".test_compute_pca_obs_sup", .test_compute_pca_obs_sup())


