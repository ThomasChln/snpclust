.actual_snp_indices <- function(object, ...) {
	dots <- list(...)
  actual_idx <- if (length(dots) == 0) {
      object@snps_idx
    } else  {
      object@snps_idx[dots[[1]]]
    }

  actual_idx
}

.actual_scan_indices <- function(object, ...) {
	dots <- list(...)
  actual_idx <- if (length(dots) == 0) {
      object@scans_idx
    } else  {
      object@scans_idx[dots[[1]]]
    }

  actual_idx
}

#' GenotypeDataSubset wrapper
#' @inheritParams GWASTools::getSnpID
#' @export
setMethod("getSnpID",
  signature(object = "GenotypeDataSubset"),
  function(object, ...) {
    callNextMethod(object, .actual_snp_indices(object, ...))
})


#' GenotypeDataSubset wrapper
#' @inheritParams GWASTools::getPosition
#' @export
setMethod("getPosition",
  signature(object = "GenotypeDataSubset"),
  function(object, ...) {
    callNextMethod(object, .actual_snp_indices(object, ...))
  })

#' GenotypeDataSubset wrapper
#' @inheritParams GWASTools::getChromosome
#' @export
setMethod("getChromosome",
  signature(object = "GenotypeDataSubset"),
  function(object, ...) {
    callNextMethod(object, .actual_snp_indices(object, ...))
  })

#' GenotypeDataSubset wrapper
#' @inheritParams GWASTools::getAlleleA
#' @export
setMethod("getAlleleA",
  signature(object = "GenotypeDataSubset"),
  function(object, ...) {
    callNextMethod(object, .actual_snp_indices(object, ...))
  })

#' GenotypeDataSubset wrapper
#' @inheritParams GWASTools::getAlleleB
#' @export
setMethod("getAlleleB",
  signature(object = "GenotypeDataSubset"),
  function(object, ...) {
    callNextMethod(object, .actual_snp_indices(object, ...))
  })

#' GenotypeDataSubset wrapper
#' @inheritParams GWASTools::getScanID
#' @export
setMethod("getScanID",
  signature(object = "GenotypeDataSubset"),
  function(object, ...) {
    callNextMethod(object, .actual_scan_indices(object, ...))
  })

#' GenotypeDataSubset wrapper
#' @inheritParams GWASTools::nsnp
#' @export
setMethod("nsnp", "GenotypeDataSubset",
  function(object) {
    length(object@snps_idx)
  })

#' GenotypeDataSubset wrapper
#' @inheritParams GWASTools::nscan
#' @export
setMethod("nscan", "GenotypeDataSubset",
  function(object) {
    length(object@scans_idx)
  })

#' GenotypeDataSubset wrapper
#' @inheritParams GWASTools::getGenotype
#' @param snp snp
#' @param scan scan
#' @export
setMethod("getGenotype",
  signature(object = "GenotypeDataSubset"),
  function(object, snp = NULL, scan = NULL, ...) {

    snps_idx <- if (is.null(snp)) {
        object@snps_idx
      } else {
        len <- if (snp[2] == -1) nsnp(object) - snp[1] + 1L else snp[2]
        object@snps_idx[seq.int(snp[1], length.out = len)]
      }

    scans_idx <- if (is.null(scan)) {
        object@scans_idx
      } else {
        len <- if (scan[2] == -1) nscan(object) - scan[1] + 1L else scan[2]
        object@scans_idx[seq.int(scan[1], length.out = len)]
      }

    res <- fetch_genotypes(as(object, 'GenotypeData'), snps_idx, scans_idx, ...)

    # bloody drop
    d <- dim(res)
    if (!is.null(d) && (d[1] == 1 || d[2] == 1)) {
      res <- as.vector(res)
    }

    res
  })


#' @export
fetch_allele1.GenotypeDataSubset <- function(obj, snps_idx) {
  idx <- if (missing(snps_idx)) obj@snps_idx else obj@snps_idx[snps_idx]

  # Force dispatch to GenotypeData for all subsequent methods
  class(obj) <- "GenotypeData"

  fetch_allele1(obj, idx)
}

#' @export
fetch_allele2.GenotypeDataSubset <- function(obj, snps_idx) {
  idx <- if (missing(snps_idx)) obj@snps_idx else obj@snps_idx[snps_idx]

  # Force dispatch to GenotypeData for all subsequent methods
  class(obj) <- "GenotypeData"

  fetch_allele2(obj, idx)
}
