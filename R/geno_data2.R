
fetch_hgdp <- function(paths) {
  urls <- c('http://www.hagsc.org/hgdp/data/hgdp.zip', 
    'ftp://ftp.cephb.fr/hgdp_v3/hgdp-ceph-unrelated.out')
  for (i in 1:2) if (!file.exists(paths[i])) download.file(urls[i], paths[i])
}

fetch_reduced_hgdp <- function(paths, zippaths, nsnp = 2e3, nscan = -1, only_eu = TRUE) {
  fetch_hgdp(paths)
  if (nsnp == -1 && nscan == -1) return()
  paths <- file.path(getwd(), paths)

  setup_temp_dir()
  txts_paths <- unzip(paths[1], zippaths, junkpaths = TRUE)
  file_paths <- c(paths[2], txts_paths[2])

# fread is bugged (since 1.9: it can not parse header with (first) empty var)
  cols <- strsplit(readLines(txts_paths[1], 1), '\t')[[1]]
  cols[1] <- 'V1'
  geno <- as.data.frame(data.table::fread(txts_paths[1], nrows = nsnp))
  names(geno) <- cols

  l_files <- list(geno)
  l_files[2:3] <- lapply(file_paths, function (i) as.data.frame(data.table::fread(i)))

  if (only_eu) l_files[[2]] <- l_files[[2]][l_files[[2]]$Region == 'Europe', ]

  if (nscan != -1) l_files[[2]] <- l_files[[2]][seq_len(nscan), ]
  idxs <- match(l_files[[2]][[1]], colnames(l_files[[1]]))
  l_files[[2]] <- l_files[[2]][!is.na(idxs), ]
  write.table(l_files[[2]], paths[2], row.names = FALSE,
    sep = '\t', quote = FALSE)
  l_files[[3]] <- l_files[[3]][match(l_files[[1]][[1]], l_files[[3]][[1]]), ]
  write.table(l_files[[3]], txts_paths[2], row.names = FALSE, col.names = FALSE,
    sep = '\t', quote = FALSE)
  l_files[[1]] <- l_files[[1]][c(1, na.omit(idxs))]
  write.table(l_files[[1]], txts_paths[1], row.names = FALSE,
    sep = '\t', quote = FALSE)

  dir.create('hgdp')
  for (i in 1:2) file.rename(txts_paths[i], zippaths[i]) 
  zip(paths[1], zippaths)
}

#' save_hgdp_as_gds
#' @param paths    Paths of the zip and txt files
#' @param zippaths Paths of the genotype and snp files in the zip
#' @return Path of the saved gds file
#' @export
save_hgdp_as_gds <- function(
  paths = file.path(system.file('extdata', package = 'snpclust'),
    paste0('hgdp.', c('zip', 'txt', 'gds'))),
  zippaths = paste0('hgdp/', c('HGDP_FinalReport_Forward.txt', 'HGDP_Map.txt'))
  ) {
  setup_temp_dir()
  txts_paths <- unzip(paths[1], zippaths, junkpaths = TRUE)
  actg_gdata <- actg_tsv_to_gdata(txts_paths[1], paths[2],
    c('scan_id', 'gender', 'population', 'geographic_origin', 'region'),
    txts_paths[2])
  file.remove(txts_paths)
  save_genotype_data_as_gds(actg_gdata, paths[3], quiet = TRUE)
  paths[3]
}

###############################################################################
gds_to_bedtargz <- function(gds, tarpath = gsub('gds$', 'tar.gz', gds)) {
  # get plink files with gds
  tarname <- gsub('[.]tar[.]gz$', '', tarpath)
  .gds_to_bedtargz(gds, tarname)

  # replace snp probe ids by ones from gdata 
  gdata <- load_gds_as_genotype_data(gds)
  on.exit(close(gdata))
  bim <- paste0(tarname, '.bim')
  df_snps <- read.table(bim, sep = '\t')
  df_snps[[2]] <- gdata@snpAnnot@data$probe_id 
  write.table(df_snps, bim, sep = '\t',
    row.names = FALSE, col.names = FALSE, quote = FALSE)

  plink_files <- paste0(tarname, '.', c('bed', 'bim', 'fam'))
  tarcmd <- paste('-czf', tarpath, paste(plink_files, collapse = ' '))
  system2('tar', tarcmd)
  file.remove(plink_files)
  tarpath
}

.gds_to_bedtargz <- function(gds, tarname) {
  gdsobj <- SNPRelate::snpgdsOpen(gds, FALSE)
  on.exit(closefn.gds(gdsobj))
  if (any(grepl('snp[.]allele', capture.output(print(gdsobj))))) {
    delete.gdsn(index.gdsn(gdsobj, 'snp.allele'))
  }
  SNPRelate::snpgdsGDS2BED(gdsobj, tarname, verbose = FALSE)
}

###############################################################################
#' actg_tsv_to_gdata
#'
#' Converts TSV files (geno, SNPs, scans) with genotypes as ACTG characters
#' to GWASTools GenotypeData.
#' Jun 4, 2014
#'
#' @param geno_path     Path of genotype text tsv file.
#'                      First line is scan ids and first column is snps ids.
#'                      NAs are specified via na_encoding parameter.
#'                      Header is managed automatically by data.table::fread.
#' @param scans_path    Path of scans text tsv file.
#'                      Header is managed automatically by data.table::fread.
#' @param scans_col_map Character vector mapping the columns of the scans file,
#'                      the index in the vector is the column index in the text
#'                      file, the characters are the names for the gds file
#'                      NAs specify undesirable columns.
#'                        Required characters: 'scan_id'
#' @param snps_path     Path of snps text tsv file
#'                      Header is managed automatically by data.table::fread.
#' @param snps_col_map  Character vector mapping the columns of the snpss file,
#'                      the index in the vector is the column index in the text
#'                      file, the characters are the names for the gds file.
#'                      NAs specify undesirable columns.
#'                        Required characters: probe_id, chromosome, position
#' @param na_encoding   Vector of characters encoding NAs in genotype file,
#' @return GenotypeData object
#'
#' @author tcharlon
actg_tsv_to_gdata <- function(geno_path, scans_path,
  scans_col_map = 'scan_id', snps_path,
  snps_col_map = c('probe_id', 'chromosome', 'position'),
  na_encoding = '--') {

  require_defaults <- function(param, param_name, fun) {
    required_snps_col <- as.character(formals(fun)[[param_name]])[-1]
    missing <- required_snps_col[!(required_snps_col %in% param)]
    if (!length(missing)) {
      paste0('Missing fields in ', param_name, ': ',
        paste0(missing, collapse = ', '))
    }
  }
  require_defaults(snps_col_map, 'snps_col_map', actg_tsv_to_gdata)
  require_defaults(scans_col_map, 'scans_col_map', actg_tsv_to_gdata)

# read all, convert snps and scans
  snps <- txt_snps_to_df(snps_path, snps_col_map)
  scans <- txt_scans_to_df(scans_path, scans_col_map)

  sep <- '\t'
  cols <- strsplit(readLines(geno_path, 1), sep)[[1]]
# can not pass by df for memory
  large_memory <- FALSE
  geno <- if (large_memory) {
# fread is bugged (since 1.9: it can not parse header with (first) empty var)
      cols[1] <- 'V1'
      geno <- data.table::fread(geno_path)
      geno <- as.data.frame(geno)
      names(geno) <- cols
      rownames(geno) <- geno[[1]]
      geno[-1]
    } else {
      geno <- scan(geno_path, 'character', sep = '\t', quiet = TRUE)
      geno <- matrix(geno, ncol = length(cols), byrow = TRUE)
      dimnames(geno) <- list(geno[, 1], geno[1, ])
      geno[-1, -1]
    }

# order genotype following SNPs and scans
  geno <- geno[match(snps$probe_id, rownames(geno)), ]
  scan_order_idxs <- match(scans$scan_id, colnames(geno))
  scans <- scans[!is.na(scan_order_idxs), ]
  scans <- cbind(scanID = seq_along(scans[[1]]), scans)
  geno <- geno[, na.omit(scan_order_idxs)]

# convert genotype and write
  geno <- t(apply(geno, 1, actg_to_numeric, na_encoding))
  stopifnot(all(geno <= 2, na.rm = TRUE))

  build_gwastools(geno, scans, snps)
}

###############################################################################
# Converts a character actg matrix to a numeric 0,1,2 matrix
actg_to_numeric <- function(snp, na_encoding = NA) {
  f_snp <- factor(snp, exclude = na_encoding, nmax = 3)
  i_snp <- as.integer(f_snp) - 1L
  if (nlevels(f_snp) == 2) {
    lvls <- levels(f_snp)
    dble_homo <- !grepl(substr(lvls[1], 1, 1), lvls[2])
    if (dble_homo) i_snp <- i_snp * 2L
  }

  i_snp
}

###############################################################################
# read and convert to factors
txt_scans_to_df <- function(path, col_map) {
  data <- as.data.frame(data.table::fread(path, '\t'))
  data <- data[!is.na(col_map)]
  names(data) <- na.omit(col_map)
  id_idx <- match('scan_id', names(data))
  data[-id_idx] <- lapply(data[-id_idx], factor)

  data
}

###############################################################################
# read, convert to integers, and order
txt_snps_to_df <- function(path, col_map) {
  data <- as.data.frame(data.table::fread(path, '\t'))
  data <- data[!is.na(col_map)]
  names(data) <- na.omit(col_map)
  data$chromosome <- as.integer(factor(data$chromosome,
      c(1:22, 'X', 'XY', 'Y', 'M'), nmax = 26, exclude = NULL))
  data$position <- as.integer(data$position)
  data <- data[order(data$chromosome, data$position), ]

  cbind(snpID = seq_along(data[[1]]), data, alleleA = NA, alleleB = NA)
}

###############################################################################
#' build_gwastools
#'
#' wrapper for GWASTools functions to build GenotypeData object
#' Aug 2, 2014
#'
#' @param geno genotype data frame
#' @param scans scans data frame
#' @param snps snps data frame
#' @return genotype data object
#'
#' @author tcharlon
build_gwastools <- function(geno, scans, snps) {
  stopifnot(inherits(geno, 'matrix'))

  rownames(scans) <- NULL
  rownames(snps) <- NULL

  if (!('snpID' %in% names(snps))) snps$snpID <- seq_len(nrow(snps))
  if (!('scanID' %in% names(scans))) scans$scanID <- seq_len(nrow(scans))
  geno <- MatrixGenotypeReader(geno, snps$snpID, snps$chromosome, snps$position,
    scans$scanID)
  snps <- SnpAnnotationDataFrame(snps)
  scans <- ScanAnnotationDataFrame(scans)

  GenotypeData(geno, snps, scans)
}

###############################################################################
#' bed_targz_to_gds
#'
#' Untar and call snpgdsBED2GDS
#' Sep 24, 2014
#'
#' @param tarpath path of input tar file
#' @param gdspath path of output gds file
#' @return NULL
#'
#' @author tcharlon
bed_targz_to_gds <- function(tarpath, gdspath) {
  check_fx_args(tarpath = '!C1', gdspath = '!C1')
  die_unless(file.exists(tarpath), paste('File does not exist:', tarpath))
  tmp_dir <- setup_temp_dir(FALSE)
  paths <- .untar_bed(tarpath, tmp_dir)

  do.call(SNPRelate::snpgdsBED2GDS,
    append(paths[c(1, 3, 2)], list(gdspath, verbose = FALSE)))
}

###############################################################################
.untar_bed <- function(tarpath, tmp_dir) {
  paths <- untar(tarpath, list = TRUE)
  patterns <- paste0('[.]', c('bed', 'bim', 'fam'), '$')
  idxs <- sapply(patterns, grep, paths)
  # if the idxs are in a list, at least one pattern is not found
  die_if(inherits(idxs, 'list'), {
      missing <- patterns[!sapply(idxs, length)]
      paste('Pattern not found in tar file:', paste(missing, collapse = ', '))
    })
  paths <- paths[idxs]
  untar(tarpath, paths, exdir = tmp_dir, compressed = TRUE)

  file.path(tmp_dir, paths)
}

###############################################################################
.gds_read_index <- function(field, gds) read.gdsn(index.gdsn(gds, field))

###############################################################################
# Format GDS sample annotations to comply with plink format
# which is expected from the function SNPRelate::snpgdsBED2GDS
# set sample.annot first two fields as sex and phenotype
.format_gds_for_BED2GDS <- function(gdspath, sample_annot, sex, phenotype) {
  arg_check <- function(value, spec) check_arg(spec, value)
  sapply(c(gdspath, sex, phenotype), arg_check, '!C1')
  gds <- snp_gds_open(gdspath, FALSE)
  on.exit(closefn.gds(gds))

  plink_cols <- paste(sample_annot, c(sex, phenotype), sep ='/')
  l_plink_info <- lapply(plink_cols, .gds_read_index, gds)
  df_plink_info <- as.data.frame(l_plink_info)
  names(df_plink_info) <- c('sex', 'phenotype')

  delete.gdsn(index.gdsn(gds, sample_annot), TRUE)
  add.gdsn(gds, 'sample.annot', df_plink_info)
}

plink_merge <- function(tar_paths, dir, outpath = 'merge') {
  if (system('plink', ignore.stdout = TRUE) == 127) stop('plink not found')
  check_fx_args(tar_paths = '!C+', outpath = '!C1')
  l_paths <- lapply(tar_paths, .untar_bed, dir)

  # get common SNPs
  snp_files <- sapply(l_paths, function(paths) grep('[.]bim$', paths))
  snps <- read.table(l_paths[[1]][snp_files[1]])[[2]]
  for (index in seq_along(l_paths)[-1]) {
    snps_match <- read.table(l_paths[[index]][snp_files[index]])[[2]]
    snps <- intersect(snps, snps_match)
  }
  write(as.character(snps), 'common_snps')

  # add dataset id to observation names
  lapply(l_paths, function(paths) {
      scan_meta <- read.table(paths[3])
      id <- gsub('.*/|[.].*', '', paths[3])
      scan_meta[[2]] <- paste0(id, '_', scan_meta[[2]])
      write.table(scan_meta, paths[3],
        quote = FALSE, row.names = FALSE, col.names = FALSE)
    })

  # write filenames and merge
  path <- gsub('[.][a-z]+$', '', l_paths[[1]][1])
  paths <- sapply(l_paths[-1], paste, collapse = ' ')
  writeLines(paste(paths, collapse = '\n'), 'file_list')
  outpath <- file.path(dir, outpath)
  cmd_merge <- paste('--bfile', path,
    '--extract common_snps --merge-list file_list --make-bed --allow-no-sex',
    '--filter-cases --out', outpath, '> /dev/null')
  system2('plink', cmd_merge)

  # allele flip if needed
  logfile <- readLines(paste0(outpath, '.log'))
  error_regexp <- 'Error: [0-9]+ variants with 3[+] alleles present'
  if (any(grepl(error_regexp, logfile))) {
    warning(logfile)
    cmd_flip <- paste0('--bfile ', path, ' --flip ', outpath,
      '-merge.missnp --make-bed --out ', path, ' > /dev/null')
    system2('plink', cmd_flip)
    system2('plink', cmd_merge)
  }

  paste0(outpath, '.', c('bed', 'bim', 'fam'))
}
