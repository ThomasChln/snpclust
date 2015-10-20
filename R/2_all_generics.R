NULL


is_snp_first_dim <- function(obj, ...) {
  UseMethod('is_snp_first_dim')
}

get_scan_annot <- function(obj, ...) {
  UseMethod('get_scan_annot')
}

get_snp_annot <- function(obj, ...) {
  UseMethod('get_snp_annot')
}



fetch_gds <- function(obj, ...) {
  UseMethod('fetch_gds')
}


fetch_allele1 <- function(obj, snps_idx) {
  UseMethod('fetch_allele1')
}

fetch_allele2 <- function(obj, snps_idx) {
  UseMethod('fetch_allele2')
}
