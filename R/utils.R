.import_gds <- function(geno, snps) {
  l_geno <- strsplit(geno, '\t')

  ## scan ids
  scans <- as.data.frame(t(sapply(l_geno, '[', c(2, 5:6))))
  colnames(scans) <- c('db_id', 'sex_fam', 'phenotype')
  stopifnot(!any(duplicated(scans$db_id)))
  scans$scanID <- 1:nrow(scans)

  ## snp ids
  colnames(snps) <- c('chromosome', 'probe_id', 'morgan', 'position')
  stopifnot(!any(duplicated(snps$V2)))

  ## genotype
  df_geno <- as.data.frame(lapply(l_geno, '[', -(1:6)))
  m_geno <- apply(df_geno, 1, function(i) as.numeric(factor(i, exclude = '0 0')))
  lvls <- apply(df_geno, 1, function(i) levels(factor(i, exclude = '0 0')))

  counts <- apply(m_geno, 2, tabulate, 3)

  # Set major allele to 0
  revs <- apply(counts, 2, function(i) i[3] > i[1])
  m_geno[, revs] <- 4L - m_geno[, revs]

  # Set homozygote to 0 in binaries
  bins <- which(!apply(counts, 2, `[`, 3))
  first_is_homo <- sapply(lvls[bins], function(i) do.call(all.equal, as.list(strsplit(i[1], ' ')[[1]])))
  revs <- bins[which(first_is_homo != 'TRUE')]
  m_geno[, revs] <- 3L - m_geno[, revs]

  # Set indels to NA
  indels <- grep('D|I', df_geno[[1]])
  m_geno[, indels] <- NA

  # Get alleles
  lvls[revs] <- lapply(lvls[revs], rev)
  alleles <- unlist(strsplit(sapply(lvls, '[', 2), ' '))
  df_alleles <- as.data.frame(matrix(alleles, ncol = 2, byrow = TRUE))
  df_alleles[indels, ] <- NA
  df_alleles[1:2] <- lapply(df_alleles[1:2], factor)
  colnames(df_alleles) <- paste0('allele', c('A', 'B'))
  snps <- cbind(snps, df_alleles)

  # Coerce genotype to integer
  dims <- dim(m_geno)
  m_geno <- as.integer(m_geno - 1)
  dim(m_geno) <- dims 

  list(geno = m_geno, scans = scans, snps = snps)
}

import_gds <- function(dirpath, gdsname) {
  path <- grep('bed', list.files(dirpath, full.names = TRUE), value = TRUE)
  system(paste('plink --bfile', gsub('.bed$', '', path),
      '--recode --tab --out pedfile'))
  on.exit(file.remove(paste0('pedfile.', c('ped', 'map', 'log'))))
  geno <- readLines('pedfile.ped')
  snps <- utils::read.table('pedfile.map')

  l_gds <- .import_gds(geno, snps)

  gdata <- build_gwastools(t(l_gds$geno), l_gds$scans, l_gds$snps)

  if (!dir.exists('data/gds')) dir.create('data/gds')
  save_genotype_data_as_gds(gdata, gdsname)
}

impute_bed <- function(l_paths, df_snp, l_peaks, scan_ids, out,
  impute = sample_impute) {
  if (length(l_paths) > 1) plink_merge(l_paths, '.', out) else {
    l_paths[[1]] %>% file.rename(
      paste0(out, gsub('.*([.][a-z]+)$', '\\1', .)))
  }
  df_snp$probe_id[sort(unique(unlist(l_peaks)))] %>%
    plink_extract(out, ., scan_ids) %>% snprelate_impute(out, impute)
}

snprelate_impute <- function(params, path, impute) {
  tmpfile = tempfile()
  do.call(SNPRelate::snpgdsBED2GDS, c(params, tmpfile, verbose = FALSE))
  gds <- openfn.gds(tmpfile, FALSE)
  geno <- .gds_read_index('genotype', gds)
  geno[geno == 3] <- NA
  write.gdsn(index.gdsn(gds, 'genotype'), impute(geno))
  closefn.gds(gds)
  SNPRelate::snpgdsGDS2BED(tmpfile, path, verbose = FALSE)
}

plink_extract <- function(path, probe_ids, scan_ids) {
  df_scans <- utils::read.table(paste0(path, '.fam'))[as.numeric(scan_ids), ]
  utils::write.table(df_scans, 'scans', sep = ' ',
    quote = FALSE, row.names = FALSE, col.names = FALSE)
  .system2('plink', paste('--bfile', path, '--snps', paste(probe_ids, collapse = ','),
      '--keep scans --make-bed --out'))
  as.list(paste0('output_plink.', c('bed', 'fam', 'bim')))
}

.system2 <- function(exe, cmd) {
  output <- paste0('output_', exe)
  system2(exe, paste(cmd, output, '| tail -n 1 > stdout'))
  stdout <- readLines('stdout')
  if (grepl('error', stdout, TRUE)) stop(stdout)
  if (!any(grepl(output, list.files()))) stop('Output file does not exist')
}

ids_subset <- function(positions, chromosome, df_snp) {
  subset <- df_snp$chromosome == chromosome &
    df_snp$position > positions[1] &
    df_snp$position < positions[2]

  df_snp$snpID[subset]
}

setup_temp_dir <- function(chdir = TRUE, ...) {
  dir <- tempfile(...)
  dir.create(dir, recursive = TRUE)
  old_dir <- NULL
  if (chdir) old_dir <- setwd(dir)

  # on one line because it not seen by the coverage
  cleanup <- bquote({if (.(chdir)) setwd(.(old_dir));unlink(.(dir), recursive = TRUE)})

  do.call(add_on_exit, list(cleanup, parent.frame()))

  invisible(normalizePath(dir))
}


.STRING <- 'string'
check_fx_args <- function(...) {
  dots <- list(...)
  formals <- formals(sys.function(sys.parent()))

  die_if(!length(dots), 'check args specifications must not be empty')
  die_unless(all(names(dots) %in% names(formals)),
    'Bad check arg specifications, unknown variables')

  env <- parent.frame()
  call <- deparse(sys.call(sys.parent()))
  call <- paste0(call, collapse = '')

  for (arg_name in names(dots)) {
    check_arg(dots[[arg_name]], get(arg_name, env), arg_name, call)
  }

  invisible()
}

check_arg <- function(spec, value, argname = deparse(substitute(value)),
  call = '') {
  if (!is.character(spec) || length(spec) != 1 || !nzchar(spec)) {
    stop('Bad spec argument:', spec)
  }

  if (!is.character(argname) || length(argname) != 1 || !nzchar(argname)) {
    stop('Bad argname argument:', argname)
  }

  if (!is.character(call) || length(call) != 1) {
    stop('Bad call:', call)
  }

  arg_fmt <- .parse_arg_spec(spec)

  # arg msg error prefix
  err <- sprintf('Error checking arg [%s] in call %s: ', argname, call)
  .msg <- function(...) {
    paste0(err, sprintf(...))
  }

  ### check null
  if (is.null(value)) {
    die_unless(arg_fmt$null, .msg('NULL not allowed'))
    return(invisible(TRUE))
  }

  ### check length
  len <- arg_fmt$length
  vl <- length(value)
  if (is.integer(len)) {
    die_unless(vl %in% len,  .msg('expected length %i, got %i', len, vl))
  } else {
    if (len == '+') {
      die_if(vl == 0, .msg('length must be > 0'))
    } else if (grepl('\\+', len)) {
      min_le <- as.integer(gsub('\\+$', '', len))
      die_if(vl < min_le, .msg('length must be > %s', min_le - 1))
    } else {

    }
  }

  ### check type
  type <- arg_fmt$type
  if (type != 'any') {
    rtype <- if (type == .STRING) 'character' else type
    die_unless(methods::is(value, rtype), .msg('expected type %s, got %s', rtype,
        class(value)))
    if (type == .STRING) {
      die_unless(all(nzchar(value)), .msg('empty strings are not allowed'))
    }
  }

  ### check NA
  if (!arg_fmt$na) { # NAs not allowed
    die_if(any(is.na(value)), .msg('no NAs are allowed'))
  }

  invisible(TRUE)
}


.parse_arg_spec <- function(fmt) {
  die_if(!nzchar(fmt), 'Bad spec: can not be ""')
  die_unless(nchar(fmt) >= 2, 'Bad spec %s', fmt)

  # parse optional !
  null <- TRUE
  if (substr(fmt, 1, 1) == '!') {
    null <- FALSE
    fmt <- substr(fmt, 2, nchar(fmt))
  }

  # parse type
  t <- substr(fmt, 1, 1)
  TYPES <- c(i = 'integer',
             n = 'numeric',
             c = 'character',
             s = .STRING,
             b = 'logical',
             l = 'list',
             d = 'data.frame',
             a = 'any')
  type <- unname(TYPES[tolower(t)])
  die_if(is.na(type), "Bad spec type identifier: %s", t)

  na <- t == tolower(t)
  l <- nchar(fmt)
  last <- substr(fmt, l, l)

#  na <- last != '!'
  len_expr <- substr(fmt, 2, l)
  if (len_expr == '*' || len_expr == '+' || grepl('^\\d+\\+$', len_expr)) {
    len <- len_expr
  } else if (len_expr == '?') {
    len <- 0:1
  } else {
    die_if(grepl('\\+', len_expr), 'Bad spec: only allowed + or [0-9]+')
    len <- eval(parse(text = len_expr))
    nn <- try(as.integer(len), silent = TRUE)
    if (is.integer(nn)) {
      len <- sort(nn)
      die_if(len < 0, "Bad spec negative length: %i", len)
    }
  }

  die_unless(len == '*' || len == '+' || is.integer(len) ||
      grepl('^\\d+\\+$', len_expr),
    'Bad spec length: %s', len)

  res <- list(type = type, length = len, na = na, null = null)

  res
}


die_unless <- function(cond, format, ...) {
  if (!is.logical(cond) || !length(cond)) {
    stop("Bad logical argument cond: ", cond)
  }
  if (missing(format)) {
    stop("Argument 'format' is missing in die_unless()")
  }

  if (any(!cond)) {
    msg <- .build_error_message(format, ...)
    stop(msg, call. = FALSE)
  }

  invisible()
}

# build the error message, using a sprintf format string and additional
# arguments. If some of the arguments are vectors, they are collapsed as
# a character string using ','
.build_error_message <- function(format, ...) {
  dots <- list(...)
  if (!length(dots)) return(format)

  .collapse_args <- function(x) {
    if (length(x) > 1) paste0(x, collapse = ',') else x
  }
  fixed_dots <- lapply(dots, .collapse_args)

  do.call(sprintf, c(list(format), fixed_dots))
}



die_if  <- function(cond, format, ...) die_unless(!cond, format, ...)


df_columns_exist <- function(data_frame, cols, silent = FALSE) {

  check_fx_args(cols = '!S+', data_frame = '!d+')

  res <- all(cols %in% colnames(data_frame))

  if (!silent && !res) {
    cols_not_in_df <- setdiff(sort(cols), sort(colnames(data_frame)))
    warning(sprintf('Column(s) "%s" not in dataframe',
        paste(cols_not_in_df, collapse = ', ')))
  }

  res
}


catch_warnings <- function(expr) {
  ws <- list()
  res <- suppressWarnings(withCallingHandlers(expr,
      warning = function(w) { ws[[length(ws)+1]] <<- w}))
  return(list(result = res, warnings = ws))
}


"%//%" <- function(a, b) {
  if (is.null(a)) b else a
}

add_on_exit <- function(expr, where = parent.frame()) {
  do.call("on.exit", list(substitute(expr), add = TRUE), envir = where)
}


