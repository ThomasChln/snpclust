

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

  # fread can not parse header with first empty var (v1.9)
  cols <- strsplit(readLines(txts_paths[1], 1), '\t')[[1]]
  cols[1] <- 'V1'
  geno <- data.table::fread(txts_paths[1], nrows = nsnp, data.table = FALSE)
  names(geno) <- cols

  l_files <- list(geno)
  l_files[2:3] <- lapply(file_paths,
    function (i) data.table::fread(i, data.table = FALSE))

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
