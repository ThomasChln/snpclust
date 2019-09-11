genotype_data_subset <- function(gdata, snps_idxs, scans_idxs) {
  qced_geno <- fetch_genotypes(gdata)[snps_idxs, scans_idxs]
  qced_scans <- gdata@scanAnnot@data[scans_idxs, ]
  qced_snps <- gdata@snpAnnot@data[snps_idxs, ]
  qced_gdata <- build_gwastools(qced_geno, qced_scans, qced_snps)
}


fetch_genotypes <- function(
  gdata,
  snps_idx = NULL,
  scans_idx = NULL,
  char = FALSE,
  snps_first = is_snp_first_dim(gdata)) {

  stopifnot(methods::is(gdata, 'GenotypeData'))

  genos <- NULL
  if (is.null(snps_idx) && is.null(scans_idx)) {
    genos <- getGenotype(gdata)
    # getGenotype may drop
    if (is.null(dim(genos))) {
      dim(genos) <- c(nsnp(gdata), nscan(gdata))
    }
  } else {
    if (is.null(snps_idx)) snps_idx <- seq_len(nsnp(gdata))
    if (is.null(scans_idx)) scans_idx <- seq_len(nscan(gdata))
    genos <- fetch_genotypes_by_idx(gdata, snps_idx, scans_idx, snps_first)
  }

  if (char) {
    if (is.null(snps_idx)) snps_idx <- seq_len(nsnp(gdata))
    if (is.null(scans_idx)) scans_idx <- seq_len(nscan(gdata))
    a1 <- fetch_allele1(gdata, snps_idx)
    a2 <- fetch_allele2(gdata, snps_idx)
    genos <- genotypeToCharacter(genos, a1, a2)
  }

  genos
}

.is_block <- function(v) {
  (utils::tail(v, 1) == v[1] + length(v) - 1L) && !is.unsorted(v)
}


.split_sorted_ints_by_blocks <- function(ints) {
  # try to be smart by testing obvious use cases first
  if (utils::tail(ints, 1) == ints[1] + length(ints) - 1) {
    return(list(c(ints[1], length(ints))))
  }

  grps <- cumsum(c(0, diff(ints) > 1))
  blks <- split(ints, grps)
  lapply(blks, function(x) c(x[1], length(x)))
}

fetch_genotypes_by_idx <- function(
  gdata,
  snps_idx = seq_len(nsnp(gdata)),
  scans_idx = seq_len(nscan(gdata)),
  snps_first = is_snp_first_dim(gdata)
) {

  snps_order <- NULL
  is_snps_block <- .is_block(snps_idx)
  snps_blocks <- if (is_snps_block) { # be smart; is it already a block ?
      list(c(snps_idx[1], length(snps_idx)))
      } else {
          snps_order <- order(snps_idx)
          .split_sorted_ints_by_blocks(snps_idx[snps_order])
      }

  scans_order <- NULL
  is_scans_block <- .is_block(scans_idx)
  scans_blocks <- if (is_scans_block) { # be smart; is it already a block ?
          list(c(scans_idx[1], length(scans_idx)))
      } else {
          scans_order <- order(scans_idx)
          .split_sorted_ints_by_blocks(scans_idx[scans_order])
      }

  genos <- fetch_genotypes_by_blocks(gdata, snps_blocks, scans_blocks, snps_first)

  if (!is_snps_block) { # need snp reordering
     if (!is_scans_block) { # need scan reordering
       genos <- genos[order(snps_order), order(scans_order), drop = FALSE]
     } else {
       genos <- genos[order(snps_order), , drop = FALSE]
     }
  } else {
    if (!is_scans_block) { # need scan reordering
      genos <- genos[, order(scans_order), drop = FALSE]
    }
  }

  genos
}


fetch_genotypes_by_blocks <- function(
  gdata,
  snps_blocks = list(c(1, nsnp(gdata))),
  scans_blocks = list(c(1, nscan(gdata))),
  snps_first = is_snp_first_dim(gdata)) {

  ### special case: one snp block - one scan block
  if (length(snps_blocks) == 1 && length(scans_blocks) == 1) {
    snp = snps_blocks[[1]]
    scan = scans_blocks[[1]]
    m <- getGenotype(gdata, snp = snp, scan = scan)
    if (is.null(dim(m))) dim(m) <- c(snp[2], scan[2])
    return(m)
  }

  nb_snps <- sum(sapply(snps_blocks, '[[', 2))
  nb_scans <- sum(sapply(scans_blocks, '[[', 2))

  # preallocate results matrix
  genotypes <- matrix(0L, nb_snps, nb_scans)
  if (snps_first) {
    col <- 1L
    for (scan in scans_blocks) {
      y <- seq.int(from = col, length.out = scan[2])
      col <- col + scan[2]
      row <- 1L
      for (snp in snps_blocks) {
        m <- getGenotype(gdata, snp = snp, scan = scan)
        if (is.null(dim(m))) dim(m) <- c(snp[2], scan[2])
        x <- seq.int(from = row, length.out = snp[2])
        row <- row + snp[2]
        genotypes[x, y] <- m
      }
    }
  } else {
    row <- 1L
    for (snp in snps_blocks) {
      x <- seq.int(from = row, length.out = snp[2])
      row <- row + snp[2]
      col <- 1L
      for (scan in scans_blocks) {
        m <- getGenotype(gdata, snp = snp, scan = scan)
        if (is.null(dim(m))) {
          dim(m) <- c(snp[2], scan[2])
        }
        y <- seq.int(from = col, length.out = scan[2])
        col <- col + scan[2]
        genotypes[x, y] <- m
      }
    }
  }

  genotypes
}

save_genotype_data_as_plink <- function(
  gdata,
  basename,
  dir = getwd(),
  snps_idx = NULL,
  scans_idx = NULL,
  ...
) {
  gdsinfo <- request_snpgds_file(gdata, snps_idx, scans_idx)
  gds <- gdsinfo$snpgds
  if (gdsinfo$new_file) {
    on.exit({closefn.gds(gds); unlink(gds$filename)}, TRUE)
  }

  SNPRelate::snpgdsGDS2BED(gds, file.path(dir, basename),
    sample.id = gdsinfo$scan_ids, snp.id = gdsinfo$snp_ids, ...)

}



save_genotype_data_as_gds <- function(
  gdata,
  filename,
  save_annotations = TRUE,
  compress = "ZIP.max",
  check = FALSE,
  quiet = FALSE,
  chunk_size = 1000,
  ...
) {
# code adapted from GWASTools::convertNcdfGds

  info <- if (quiet) function(...) {} else function(...) message(sprintf(...))

  ### snp gds

  tt <- system.time({
    allele1 <- fetch_allele1(gdata) %//% 'A'
    allele2 <- fetch_allele2(gdata) %//% 'B'

    alleles <- paste(allele1, allele2, sep = '/')
    alleles[is.na(allele1) | is.na(allele2)] <- ''
  }, gcFirst = FALSE)
  info('fetched alleles in %.2fs', tt[3])


  tt <- system.time({
    snp_infos <- data.frame(
      snpID = getSnpID(gdata),
      chromosome = getChromosome(gdata),
      position = getPosition(gdata),
      alleles = alleles
      , stringsAsFactors = FALSE)
  }, gcFirst = FALSE)
  info('fetched snp_infos in %.2fs', tt[3])

  info('saving genotypes...')
  tt <- system.time({
    it <- fetch_genotypes_by_chunk_iterator(gdata, chunk_size = chunk_size, ...)
    write_snp_gds(filename, it,
      sample_ids = getScanID(gdata), snp_infos = snp_infos, quiet = quiet)
  }, gcFirst = FALSE)
  info('took %s', tt[3])

  ### annotations
  if (save_annotations) {
    tt <- system.time({
      gds <- snp_gds_open(filename, FALSE)
      .gds_serialize_object(gds, 'snp_annot.serialized', get_snp_annot(gdata))
      .gds_serialize_object(gds, 'scan_annot.serialized', get_scan_annot(gdata))
      closefn.gds(gds)
    }, gcFirst = FALSE)
    info('saved annotations in %.2fs', tt[3])
  }

  if (check) {
    tt <- system.time({
        check_snp_gds(filename)
      }, gcFirst = FALSE)
    info('checked in %.2fs', tt[3])
  }

  invisible()
}


.gds_serialize_object <- function(gds, name, object, compress = "ZIP.max") {
  ser <- serialize(object, NULL)
  add.gdsn(gds, name, as.integer(ser), storage = "uint8", compress = compress)
}


.gds_unserialize_object <- function(gds, name) {
  ser <- read.gdsn(index.gdsn(gds, name))
  unserialize(as.raw(ser))
}



write_genotypes_to_gds <- function(
  gds,
  geno_iterator,
  node_name = 'genotype',
  quiet = FALSE)
{
  info <- if (quiet) function(...) {} else function(...) message(sprintf(...))

  .convert_chunk <- function(x) {
    x[is.na(x)] <- 3L
    t(x)
  }

  nb_chunks <- geno_iterator$nb

  ### save the first chunk
  chunk <- .convert_chunk(geno_iterator$get(1))
  node <- add.gdsn(gds, node_name, chunk, storage = "bit2")

  current <- 0L
  for (i in seq.int(from = 2, length.out = nb_chunks - 1)) {
    perc <- ceiling(i / nb_chunks * 100)
    if (perc > current) {
      info('writing genotypes %i%%', perc)
      current <- perc
    }

    chunk <- .convert_chunk(geno_iterator$get(i))
    append.gdsn(node, chunk)
  }

  put.attr.gdsn(node, "sample.order")

  node
}


# workaround for current bug when forking in gdsfmt
.create_gds <- function(filename) {
  gds <- createfn.gds(filename)
  closefn.gds(gds)
  openfn.gds(filename, readonly = FALSE, allow.duplicate = FALSE,
    allow.fork = TRUE)
}

write_snp_gds <- function(filename, geno_iterator, sample_ids, snp_infos,
  compress = "ZIP.max", quiet = FALSE, ...) {

  stopifnot(!anyDuplicated(sample_ids))
  cols <- c('snpID', 'chromosome', 'position', 'alleles')
  stopifnot(all(cols %in% colnames(snp_infos)), !anyDuplicated(snp_infos$snpID))

  info <- if (quiet) function(...) {} else function(...) message(sprintf(...))

  gfile <- .create_gds(filename)
  on.exit(closefn.gds(gfile), TRUE)

  .add <- function(name, val) {
    info('writing %s...', name)
    add.gdsn(gfile, name, val, compress = compress, closezip = TRUE)
  }

  .add("sample.id", sample_ids)
  .add("snp.id", snp_infos$snpID)
  .add("snp.position", snp_infos$position)
  .add("snp.chromosome", snp_infos$chromosome)
  .add("snp.allele", snp_infos$alleles)

  write_genotypes_to_gds(gfile, geno_iterator, quiet = quiet, ...)
}



load_gds_as_genotype_data <- function(gds_file, read_snp_annot = TRUE,
  read_scan_annot = TRUE) {

  gds <- snp_gds_open(gds_file, TRUE, TRUE, TRUE)

  .read_node <- function(name) {
    if (is.null(index.gdsn(gds, name, silent = TRUE)))
      return(NULL)
    .gds_unserialize_object(gds, name)
  }

  snp_annot <- if (read_snp_annot)
      .read_node('snp_annot.serialized')
    else
      NULL

  scan_annot <- if (read_scan_annot)
      .read_node('scan_annot.serialized')
    else
      NULL

  reader <- gds_genotype_reader(gds)
  GenotypeData(reader, snp_annot, scan_annot)
}


# temp wrapper for snpgdsOpen while it is not mainstream
snp_gds_open <- function(...) {
  SNPRelate <- NULL; Library <- library; Library(SNPRelate)
  gds <- SNPRelate::snpgdsOpen(...)

  # hack: some GWASTools functions check for gds.class
  class(gds) <- "gds.class"

  gds
}



request_snpgds_file <- function(
    gdata,
    snps_idx = NULL,
    scans_idx = NULL
) {
  gds <- fetch_gds(gdata)
  new_file <- NULL
  if (is.null(gds)) {
    # no gds, must create one
    gdsfn <- tempfile(fileext = '.gds')

    # subset if needed
    gdata2 <- if (is.null(snps_idx) && is.null(scans_idx)) {
        gdata
      } else {
        if (is.null(snps_idx)) snps_idx = seq_len(nsnp(gdata))
        if (is.null(scans_idx)) scans_idx = seq_len(nscan(gdata))
        genotype_data_subset(gdata, snps_idx, scans_idx)
      }

    save_genotype_data_as_gds(gdata2, gdsfn, quiet = TRUE)
    gds <- snp_gds_open(gdsfn, TRUE, TRUE, TRUE)
    new_file <- TRUE
    # gdata already subsetted now
    snps_idx <- scans_idx <- NULL
    snp_ids <- scan_ids <- NULL

  } else {
    gdsfn <- gds$filename
    new_file <- FALSE
    snp_ids <- if (is.null(snps_idx)) NULL else {
              getSnpID(gdata, snps_idx)
    }
    scan_ids <- if (is.null(scans_idx)) NULL else {
              getScanID(gdata, scans_idx)
          }

  }

  list(snpgds = gds, snps_idx = snps_idx, scans_idx = scans_idx,
      new_file = new_file, snp_ids = snp_ids, scan_ids = scan_ids)
}


gds_genotype_reader <- function(gds, ...) {
  if (is.character(gds)) {
    gds <- snp_gds_open(gds, TRUE, TRUE, TRUE)
  }

  genotypeDim <- if (is_snp_first_dim(gds)) "snp,scan" else "scan,snp"
  gdsreader <- methods::new('GdsReader', filename = gds$filename, handler = gds)

  methods::new('GdsGenotypeReader', gdsreader, genotypeDim = genotypeDim, ...)
}


fetch_genotypes_by_chunk_iterator <- function(
  gdata,
  chunk_size = 100L,
  nb_chunks = NULL) {

  nb_snps <- nsnp(gdata)

  if (is.null(nb_chunks)) {
    nb_chunks <- as.integer(ceiling(nb_snps / chunk_size))
  } else {
    chunk_size <- as.integer(ceiling(nb_snps / nb_chunks))
  }

  .get_chunk <- function(n) {
    die_unless(n > 0 && n <= nb_chunks, 'bad chunk index %i', n)
    i <- (n - 1) * chunk_size + 1L
    j <- min(nb_snps, i + chunk_size - 1L)
    idx <- seq.int(i, j)
    fetch_genotypes_by_idx(gdata, idx)
  }

  list(nb = nb_chunks, chunk_size = chunk_size, get = .get_chunk)
}


check_snp_gds <- function(gds) {
  if (is.character(gds)) {
    gds <- snp_gds_open(gds, TRUE, TRUE, TRUE)
    on.exit(closefn.gds(gds))
  }

  nb_snps <- objdesp.gdsn(index.gdsn(gds, "snp.id"))$dim
  die_unless(nb_snps > 0, 'bad number of snps: %s', nb_snps)
  nb_scans <- objdesp.gdsn(index.gdsn(gds, "sample.id"))$dim
  die_unless(nb_scans > 0, 'bad number of scans: %s', nb_scans)

  snp_first_dim <- is_snp_first_dim(gds)
  die_unless(is.finite(snp_first_dim), 'no genotype order information')

  geno_dim <- objdesp.gdsn(index.gdsn(gds, "genotype"))$dim
  if (!snp_first_dim) geno_dim <- rev(geno_dim)
  die_unless(identical(geno_dim, c(nb_snps, nb_scans)),
    'genotype matrix dimension mismatch')

  snp_ids <- read.gdsn(index.gdsn(gds, "snp.id"))
  die_if(anyDuplicated(snp_ids), 'snp ids not unique')

  sample_ids <- read.gdsn(index.gdsn(gds, "sample.id"))
  die_if(anyDuplicated(sample_ids), 'sample ids not unique')

  # chromosome
  chr <- read.gdsn(index.gdsn(gds, "snp.chromosome"))
  die_unless(length(chr) == nb_snps, 'bad chromosomes length')
  if (!all(is.na(chr))) {
    r <- range(chr, na.rm = TRUE)
    die_if(min(r) < 1L, 'bad chromosomes, found < 1')
    die_if(max(r) > 27L, 'bad chromosomes, found > 27')
    die_if(is.unsorted(chr, na.rm = TRUE), 'chromosomes not sorted')
  }

  # positions
  pos <- read.gdsn(index.gdsn(gds, "snp.position"))
  die_unless(length(pos) == nb_snps, 'bad positions length')
  pmin <- min(pos, na.rm = TRUE)
  die_if(pmin < 1L, 'bad positions, found < 1')
  by_chr_snp <- order(chr, pos)
  die_if(is.unsorted(by_chr_snp), 'not sorted by chromosome/position')

  # alleles
  alleles <- read_snp_gds_alleles(gds)
  if (!is.null(alleles)) {
    alleles <- unique(alleles)
    alleles <- alleles[!is.na(alleles)]
    alls <- unique(unlist(strsplit(unique(alleles), "/")))
    good_alleles <- c('A', 'C', 'T', 'G', '-', NA)
    die_unless(all(alls %in% good_alleles), 'bad alleles')
  }
}

read_snp_gds_alleles <- function(gds) {
  if (is.null(index.gdsn(gds, "snp.allele", silent = TRUE)))
    return(NULL)

  alleles <- read.gdsn(index.gdsn(gds, "snp.allele"))
  alleles[!nzchar(alleles)] <- NA

  alleles
}


check_genotype_data <- function(gdata) {

  nb_snps <- nsnp(gdata)
  die_unless(nb_snps > 0, 'gdata: bad number of snps: %s', nb_snps)
  nb_scans <- nscan(gdata)
  die_unless(nb_scans > 0, 'gdata: bad number of scans: %s', nb_scans)
  die_if(anyDuplicated(getSnpID(gdata)), 'gdata: snpIDs not unique')
  die_if(anyDuplicated(getScanID(gdata)), 'gdata: scanIDs not unique')

  ### the genotype reader
  reader <- gdata@data
  die_unless(nb_snps == nsnp(reader), 'reader: bad number of snps')
  die_unless(nb_scans == nscan(reader), 'reader: bad number of scans')

  .check_chromosome(reader, 'reader')
  die_unless(identical(getSnpID(reader), getSnpID(gdata)), 'reader: bad snpIDs')
  die_unless(identical(getScanID(reader), getScanID(gdata)), 'reader: bad scanIDs')
  die_unless(identical(getPosition(reader), getPosition(gdata)),
    'reader: bad getPosition')
  die_unless(identical(getChromosome(reader), getChromosome(gdata)),
    'reader: bad getChromosome')
  die_unless(identical(getScanID(reader), getScanID(gdata)), 'reader: bad scanIDs')

  # chromosome
  chr <- getChromosome(reader)
  chr <- chr[!is.na(chr)]
  if (length(chr) > 0) {
    r <- range(chr)
    die_if(min(r) < 1L, 'bad chromosomes, found < 1')
    die_if(max(r) > 27L, 'bad chromosomes, found > 27')
    die_if(is.unsorted(chr), 'chromosomes not sorted')
  }
  # positions
  pos <- getPosition(reader)
  pos <- pos[!is.na(pos)]
  if (length(pos) > 0) {
    die_if(min(pos) < 1L, 'bad positions, found < 1')
  }

  if (length(chr) > 0 && length(pos) > 0) {
      by_chr_snp <- order(getChromosome(reader), getPosition(reader))
      die_if(is.unsorted(by_chr_snp), 'not sorted by chromosome/position')
  }

  ### consistency with snpAnnot
  if (!is.null(sa <- gdata@snpAnnot)) {
    die_unless(identical(getSnpID(sa), getSnpID(gdata)), 'snpAnnot: bad snpIDs')
    .check_chromosome(sa, 'snpAnnot')
    die_unless(identical(getPosition(sa), getPosition(gdata)),
      'snpAnnot: bad getPosition')
    die_unless(identical(getChromosome(sa), getChromosome(gdata)),
      'snpAnnot: bad getChromosome')

    # alleles
    if (hasVariable(sa, sa@alleleACol)) {
      aa <- getAlleleA(sa)
      ab <- getAlleleB(sa)
      good_alleles <- c('A', 'C', 'T', 'G', '-', NA)
      die_unless(all(unique(aa) %in% good_alleles), 'bad alleleA')
      die_unless(all(unique(aa) %in% good_alleles), 'bad alleleB')
    }
  }

  ### consistency with scanAnnot
  if (!is.null(sa <- gdata@scanAnnot)) {
    die_unless(identical(getScanID(sa), getScanID(gdata)), 'scanAnnot: bad scanIDs')
    if (hasSex(sa)) {
      sex <- getSex(sa)
      die_unless(all(sex %in% c("F", "M", NA)), 'scanAnnot: bad getSex')
    }
  }
}


.check_chromosome <- function(o, msg) {
  die_unless(identical(autosomeCode(o), 1:22), '%s: bad autosome codes', msg)
  die_unless(XchromCode(o) == 23L, '%s: bad X chr code', msg)
  die_unless(XYchromCode(o) == 24L, '%s: bad XY chr code', msg)
  die_unless(YchromCode(o) == 25L, '%s: bad Y chr code', msg)
  die_unless(MchromCode(o) == 26L, '%s: bad Y chr code', msg)
}


