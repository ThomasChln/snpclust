
fetch_hgdp <- function(paths) {
  urls <- c('http://www.hagsc.org/hgdp/data/hgdp.zip', 
    'ftp://ftp.cephb.fr/hgdp_v3/hgdp-ceph-unrelated.out')
  for (i in 1:2) {
    if (!file.exists(paths[i])) utils::download.file(urls[i], paths[i])
  }
}

reduce_hgdp <- function(paths,
  zippaths = paste0('hgdp/', c('HGDP_FinalReport_Forward.txt', 'HGDP_Map.txt')),
  nsnp = -1, nscan = -1, region_selection = 'Europe') {

  setup_temp_dir()
  txts_paths <- utils::unzip(paths[1], zippaths, junkpaths = TRUE)
  file_paths <- c(paths[2], txts_paths[2])

# fread is bugged (since 1.9: it can not parse header with (first) empty var)
  cols <- strsplit(readLines(txts_paths[1], 1), '\t')[[1]]
  cols[1] <- 'V1'
  geno <- data.table::fread(txts_paths[1], nrows = nsnp, data.table = FALSE)
  names(geno) <- cols

  l_files <- list(geno)
  l_files[2:3] <- lapply(file_paths, function (i) data.table::fread(i, data.table = FALSE))

  if (!is.null(region_selection)) {
    l_files[[2]] <- l_files[[2]][grep(region_selection, l_files[[2]]$Region), ]
  }

  if (nscan != -1) l_files[[2]] <- l_files[[2]][seq_len(nscan), ]
  idxs <- match(l_files[[2]][[1]], colnames(l_files[[1]]))
  l_files[[2]] <- l_files[[2]][!is.na(idxs), ]
  utils::write.table(l_files[[2]], paths[2], row.names = FALSE,
    sep = '\t', quote = FALSE)
  l_files[[3]] <- l_files[[3]][match(l_files[[1]][[1]], l_files[[3]][[1]]), ]
  utils::write.table(l_files[[3]], txts_paths[2], row.names = FALSE, col.names = FALSE,
    sep = '\t', quote = FALSE)
  l_files[[1]] <- l_files[[1]][c(1, stats::na.omit(idxs))]
  utils::write.table(l_files[[1]], txts_paths[1], row.names = FALSE,
    sep = '\t', quote = FALSE)

  dir.create('hgdp')
  for (i in 1:2) file.rename(txts_paths[i], zippaths[i]) 
  utils::zip(paths[1], zippaths)
}

#' save_hgdp_as_gds
#' @param paths    Paths of the zip, txt, and gds files
#' @param zippaths Paths of the genotype and snp files in the zip
#' @return Path of the saved gds file
#' @export
save_hgdp_as_gds <- function(
  paths = file.path(system.file('extdata', package = 'snpclust'),
    paste0('hgdp.', c('zip', 'txt', 'gds'))),
  zippaths = paste0('hgdp/', c('HGDP_FinalReport_Forward.txt', 'HGDP_Map.txt'))
  ) {
  setup_temp_dir()
  txts_paths <- utils::unzip(paths[1], zippaths, junkpaths = TRUE)
  actg_gdata <- actg_tsv_to_gdata(txts_paths[1], paths[2],
    c('scan_id', 'gender', 'population', 'geographic_origin', 'region'),
    txts_paths[2])
  file.remove(txts_paths)
  save_genotype_data_as_gds(actg_gdata, paths[3], quiet = TRUE)
  paths[3]
}

gds_to_bed <- function(gds) {
  name = gsub('.gds', '', gds)
  snprelate_gds_to_bed(name)

  # replace snp probe ids by ones from gdata 
  gdata <- load_gds_as_genotype_data(gds)
  on.exit(close(gdata))
  bim <- paste0(name, '.bim')
  df_snps <- utils::read.table(bim, sep = '\t')
  df_snps[[2]] <- gdata@snpAnnot@data$probe_id 
  utils::write.table(df_snps, bim, sep = '\t',
    row.names = FALSE, col.names = FALSE, quote = FALSE)

  lapply(name, paste0, '.', c('bed', 'bim', 'fam'))
}

# get plink files with gds
snprelate_gds_to_bed <- function(name) {
  gdsobj <- paste0(name, '.gds') %>% SNPRelate::snpgdsOpen(FALSE)
  on.exit(closefn.gds(gdsobj))
  if (any(grepl('snp[.]allele', utils::capture.output(print(gdsobj))))) {
    delete.gdsn(index.gdsn(gdsobj, 'snp.allele'))
  }
  SNPRelate::snpgdsGDS2BED(gdsobj, name, verbose = FALSE)
}

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
      geno <- data.table::fread(geno_path, data.table = FALSE)
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
  geno <- geno[, stats::na.omit(scan_order_idxs)]

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
  names(data) <- stats::na.omit(col_map)
  id_idx <- match('scan_id', names(data))
  data[-id_idx] <- lapply(data[-id_idx], factor)

  data
}

###############################################################################
# read, convert to integers, and order
txt_snps_to_df <- function(path, col_map) {
  data <- as.data.frame(data.table::fread(path, '\t'))
  data <- data[!is.na(col_map)]
  names(data) <- stats::na.omit(col_map)
  data$chromosome <- as.integer(factor(data$chromosome,
      c(1:22, 'X', 'XY', 'Y', 'M'), nmax = 26, exclude = NULL))
  data$position <- as.integer(data$position)
  data <- data[order(data$chromosome, data$position), ]

  cbind(snpID = seq_along(data[[1]]), data, alleleA = NA, alleleB = NA)
}

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

bed_to_gds <- function(paths, gds) {
  do.call(SNPRelate::snpgdsBED2GDS,
    append(paths[c(1, 3, 2)], list(gdspath, verbose = FALSE)))
}

.gds_read_index <- function(field, gds) read.gdsn(index.gdsn(gds, field))

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

# get common SNPs
write_common_snps = function(l_paths) {
  snp_files <- sapply(l_paths, function(paths) grep('[.]bim$', paths))
  snps <- utils::read.table(l_paths[[1]][snp_files[1]])[[2]]
  for (index in seq_along(l_paths)[-1]) {
    snps_match <- utils::read.table(l_paths[[index]][snp_files[index]])[[2]]
    snps <- intersect(snps, snps_match)
  }
  write(as.character(snps), 'common_snps')
}

# add dataset id to observation names
write_scan_dataset_id = function(paths) {
  scan_meta <- utils::read.table(paths[3])
  id <- gsub('.*/|[.].*', '', paths[3])
  scan_meta[[2]] <- paste0(id, '_', scan_meta[[2]])
  utils::write.table(scan_meta, paths[3],
    quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# write filenames and merge
plink_merge_beds <- function(l_paths, dir, outpath) {
  path <- gsub('[.][a-z]+$', '', l_paths[[1]][1])
  paths <- sapply(l_paths[-1], paste, collapse = ' ')
  writeLines(paste(paths, collapse = '\n'), 'file_list')
  outpath <- file.path(dir, outpath)
  cmd_merge <- paste('--bfile', path,
    '--extract common_snps --merge-list file_list --make-bed --allow-no-sex',
    '--filter-cases --out', outpath, '> /dev/null')
  system2('plink', cmd_merge)
  cmd_merge
}

plink_merge <- function(l_paths, dir, outpath) {
  if (system('plink', ignore.stdout = TRUE) == 127) stop('plink not found')
  write_common_snps(l_paths)
  lapply(l_paths, write_scan_dataset_id)
  plink_merge_beds(l_paths, dir, outpath) %>% plink_allele_flip(outpath, path, .)
  paste0(outpath, '.', c('bed', 'bim', 'fam'))
}

plink_allele_flip = function(outpath, path, cmd_merge) {
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
}
