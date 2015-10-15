### maybe put a default function, or use the union type, or define it
### on MatrixGenotypeReader and NcdfGenotypeReader


## ' for GenotypeData
## ' @inheritParams get_max_snp_chunk_size
## ' @export
#get_max_snp_chunk_size.GenotypeData <- function(obj, ...) {
#  nsnp(obj)
#}
#
## ' in a GenotypeDataSubset, there are no physical chunks
## ' @inheritParams get_max_snp_chunk_size
## ' @export
#get_max_snp_chunk_size.GenotypeDataSubset <- function(obj, ...) {
#  1L
#}


fetch_gds.default <- function(obj, ...) {
  NULL
}

fetch_gds.GdsGenotypeReader <- function(obj, ...) {
  obj@handler
}

fetch_gds.GenotypeData <- function(obj, ...) {
  fetch_gds(obj@data)
}

fetch_gds.GenotypeDataSubset <- function(obj, ...) {
  NULL
}


is_snp_first_dim.default <- function(obj, ...) {
  NA
}

is_snp_first_dim.gds.class <- function(obj, ...) {
  snpfirstdim <- NA
  rd <- names(get.attr.gdsn(index.gdsn(obj, "genotype")))
	if ("snp.order" %in% rd) snpfirstdim <- TRUE
	if ("sample.order" %in% rd) snpfirstdim <- FALSE

  snpfirstdim
}


is_snp_first_dim.GenotypeData <- function(obj, ...) {
  is_snp_first_dim(obj@data)
}


is_snp_first_dim.GdsGenotypeReader <- function(obj, ...) {
  is_snp_first_dim(obj@handler)
}


is_snp_first_dim.MatrixGenotypeReader <- function(obj, ...) {
  TRUE
}

is_snp_first_dim.NcdfGenotypeReader <- function(obj, ...) {
  TRUE
}





fetch_allele1.default <- function(obj, snps_idx) { NULL }

fetch_allele2.default <- function(obj, snps_idx) { NULL }

fetch_allele1.GenotypeData <- function(obj, ...) {
  alleles <- getAlleleA(obj, ...)
  if (!is.null(alleles)) return(alleles)

  fetch_allele1(obj@data, ...)
}

fetch_allele2.GenotypeData <- function(obj, ...) {
  alleles <- getAlleleB(obj, ...)
  if (!is.null(alleles)) return(alleles)

  fetch_allele2(obj@data, ...)
}


fetch_snpgds_alleles <- function(gdsobj, snps_idx) {
  node <- index.gdsn(gdsobj, "snp.allele", silent = TRUE)
  if (is.null(node)) return(NULL)

  if (missing(snps_idx)) {
    return(read.gdsn(node))
  }

  sel <- rep(FALSE, objdesp.gdsn(node)$dim)
  sel[snps_idx] <- TRUE
  alleles <- readex.gdsn(node, sel)
  # now restore order !
  alleles[order(snps_idx)]
}

fetch_allele1.GdsGenotypeReader <- function(obj, snps_idx) {
  alleles <- fetch_snpgds_alleles(obj@handler, snps_idx)
  if (is.null(alleles)) return(NULL)
  substr(alleles, 1, 1)
}

fetch_allele2.GdsGenotypeReader <- function(obj, snps_idx) {
  alleles <- fetch_snpgds_alleles(obj@handler, snps_idx)
  if (is.null(alleles)) return(NULL)
  substr(alleles, 3, 3)
}


get_scan_annot.GenotypeData <- function(obj, ...) {
  obj@scanAnnot
}

get_snp_annot.GenotypeData <- function(obj, ...) {
  obj@snpAnnot
}

get_scan_annot.GenotypeDataSubset <- function(obj, ...) {
  obj@scanAnnot[obj@scans_idx, ]
}

get_snp_annot.GenotypeDataSubset <- function(obj, ...) {
  obj@snpAnnot[obj@snps_idx, ]
}

