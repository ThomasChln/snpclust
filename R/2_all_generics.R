#' @name All generics
#' @docType methods
#' @rdname generics
NULL


#' test if the genotypes first dimension is snp
#'
#' The genotypes matrix can either be snp*sample or sample*snp
#' It is more efficient is the matrix is to be accessed by snp to
#' store it  sample*snp, i.e. is_gds_snp_first_dim == FALSE
#'
#' @param obj   an object containing genotypes
#' @param ...   additional arguments
#' @return TRUE iff snp is the first dimension
#' @export
#' @rdname generics
is_snp_first_dim <- function(obj, ...) {
  UseMethod('is_snp_first_dim')
}

#
## ' get the maximum usable genotype chunk size
## '
## ' @param obj   an object containing genotypes
## ' @param ...   additional arguments
## ' @export
#get_max_snp_chunk_size <- function(obj, ...) {
#  UseMethod('get_max_snp_chunk_size')
#}

#' accessor for the scanAnnot slot of GenotypeData
#'
#' @param obj   an GenotypeData like object
#' @param ...   additional arguments
#' @return a ScanAnnotationDataFrame object, or NULL
#' @export
get_scan_annot <- function(obj, ...) {
  UseMethod('get_scan_annot')
}

#' accessor for the snpAnnot slot of GenotypeData
#'
#' @param obj   an GenotypeData like object
#' @param ...   additional arguments
#' @return a SnpAnnotationDataFrame object, or NULL
#' @export
get_snp_annot <- function(obj, ...) {
  UseMethod('get_snp_annot')
}



#' get the GDS object associated with a GenotypeData object if any
#'
#' @export
#' @rdname generics
fetch_gds <- function(obj, ...) {
  UseMethod('fetch_gds')
}


#' fetch the first alleles
#'
#' @param snps_idx    the indices of the snps for which to fetch
#' @return a one-character vector of alleles, or NULL if no alleles
#' @export
#' @rdname generics
fetch_allele1 <- function(obj, snps_idx) {
  UseMethod('fetch_allele1')
}

#' fetch the second alleles
#'
#' @return a one-character vector of alleles
#' @export
#' @rdname generics
fetch_allele2 <- function(obj, snps_idx) {
  UseMethod('fetch_allele2')
}