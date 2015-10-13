

###############################################################################
.untar_impute_bed <- function(tar_paths, df_snp, l_peaks, scan_ids) {
  if (length(tar_paths) > 1) {
    plink_merge(tar_paths, '.', 'bedfile')
  } else {
    paths <- .untar_bed(tar_paths, '.')
    file.rename(paths, gsub('.*[.]', 'bedfile.', paths))
  }
  probe_ids <- df_snp$probe_id[sort(unique(unlist(l_peaks)))]

  .impute_bed('bedfile', probe_ids, scan_ids)
}

.impute_bed <- function(path, probe_ids, scan_ids, impute = sample_impute) {
  df_scans <- read.table(paste0(path, '.fam'))[as.numeric(scan_ids), ]
  write.table(df_scans, 'scans', sep = ' ',
    quote = FALSE, row.names = FALSE, col.names = FALSE)
  .system2('plink', paste('--bfile', path,
      '--snps', paste(probe_ids, collapse = ','),
      '--keep scans --make-bed --out'))
  params <- as.list(paste0('output_plink.', c('bed', 'fam', 'bim')))
  do.call(SNPRelate::snpgdsBED2GDS, c(params, 'gds', verbose = FALSE))

  gds <- openfn.gds('gds', FALSE)
  geno <- .gds_read_index('genotype', gds)
  geno[geno == 3] <- NA
  write.gdsn(index.gdsn(gds, 'genotype'), impute(geno))
  closefn.gds(gds)

  SNPRelate::snpgdsGDS2BED('gds', path, verbose = FALSE)
}

.system2 <- function(exe, cmd) {
  output <- paste0('output_', exe)
  system2(exe, paste(cmd, output, '| tail -n 1 > stdout'))
  stdout <- readLines('stdout')
  if (grepl('error', stdout, TRUE)) stop(stdout)
  if (!any(grepl(output, list.files()))) stop('Output file does not exist')
}

.untar_bed_targz <- function(path) {
  untar(path, compressed = TRUE)
  paths <- untar(path, list = TRUE)
  gsub('[.].*', '', paths[grep('[.]bed$', paths)])
}

###############################################################################
ids_subset <- function(positions, chromosome, df_snp) {
  subset <- df_snp$chromosome == chromosome &
    df_snp$position > positions[1] &
    df_snp$position < positions[2]

  df_snp$snpID[subset]
}

###############################################################################
#' create a temp directory, that will automatically deleted
#'
#' when the calling function exits, and can set the working directory
#' temporarily into this new directory. The previous current directory will
#' be restored when the calling function exits.
#'
#' @param chdir   whether to set the current directory to the new temp dir
#' @param ...     args forwarded to \code{\link{tempfile}}
#' @return the created directory path invisibly
#' @export
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
#' check the function arguments using a concise syntax
#'
#' Use the \code{\link{check_arg}} syntax for each checked argument
#'
#' @param ...         specifications for arguments
#'
#' @return NULL
#' @export
#' @family args_checking
#' @examples
#' f <- function(x, y = 1, z = NULL) {
#'    check_fx_args(x = "I+", y = "n1", z = "C*")
#' }
#'
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

#' check that a value is conform to our argument format
#'
#' die if the argument is invalid
#'
#' The format string is made of
#'  \describe{
#'      \item{not_null:optional}{ a ! means that NULL is not allowed }
#'
#'      \item{type:mandatory}{ a one letter type descriptor. upper-case means
#'          that NAs are not allowed.
#'       \describe{
#'          \item{i}{integer}
#'          \item{n}{numeric}
#'          \item{c}{character}
#'          \item{s}{string: non empty strings}
#'          \item{b}{logical (boolean) }
#'          \item{l}{list}
#'          \item{d}{data.frame}
#'          \item{a}{any type}
#'        }}
#'
#'      \item{length:mandatory}{
#'          N.B:For a data.frame, the length is the number of columns
#'          does not apply to NULL value, use ! to forbid NULL
#'       \describe{
#'          \item{*}{ any length }
#'          \item{+}{ length > 0 }
#'          \item{?}{ 0 or 1 }
#'          \item{number+}{ length equal or higher than the indicated number}
#'          \item{number expression}{ a R expression that must evaluate to
#'              something convertible to integer }
#'        }}
#'
#'}
#'
#' @param spec          the check specification string
#' @param value         the value to check against the format
#' @param argname       the name of the argument to check
#' @param call          the function call, for error message clarity
#'
#' @return TRUE
#' @family args_checking
#' @export
#' @examples
#' check_arg('N+', 1:5)
#' check_arg('!S+', "coucou")
#' check_arg('i1:10', 1:5)
#' check_arg('i2+', 1:5)
#'
#'
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
    die_unless(is(value, rtype), .msg('expected type %s, got %s', rtype,
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


#' die (stop) unless the condition is TRUE, with a \code{\link{sprintf}}
#' formatted message
#'
#' @param cond     the condition that should be met
#' @param format   the formatted message
#' @param ...     additional arguments for \code{\link{sprintf}}.
#'                If some arguments have length > 1, they are automatically
#'                collapsed into a string
#'
#' @return NULL
#' @family exceptions
#' @export
#' @examples \dontrun{
#'  s <- c('A', 'b', 'C', 'd')
#'  res <- try(die_unless(s %in% LETTERS, 'bad letters: "%s"',
#'       s[!s%in% LETTERS]), silent = TRUE)
#'  # "Error : bad letters: \"b,d\"\n"
#'
#'
#' }
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



#' die (stop) if the condition is TRUE, with a \code{\link{sprintf}}
#' formatted message
#'
#' @param cond     the condition that should be met
#' @param format   the formatted message
#' @param ...     additional args for \code{\link{sprintf}}
#'
#' @return NULL
#' @family exceptions
#' @export
die_if  <- function(cond, format, ...) die_unless(!cond, format, ...)


#' Checks if one or more column names exist in the data.frame
#'
#' TRUE iff column found + warning with name(s) of missing column(s)
#'
#' @param data_frame            the data frame to check columns in
#' @param cols                  non empty vector of column names
#' @param silent                display warning (default "TRUE")
#' @return boolean
#' @export
#' @family data.frame
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


#' suppress the warnings, and return them along with the result of the expression
#'
#' @param expr 		expression to evaluate.
#'
#' @return a list with 2 elements: the result of the evaluated expression,
#'    and the list of caught warning messages
#' @family exceptions
#' @export
catch_warnings <- function(expr) {
  ws <- list()
  res <- suppressWarnings(withCallingHandlers(expr,
      warning = function(w) { ws[[length(ws)+1]] <<- w}))
  return(list(result = res, warnings = ws))
}


#' an imitation of the perl // operator
#' 
#' like "NULL OR" 
#'
#' @param   a   any object
#' @param   b   any object
#' @return  a unless NULL, in which case b
#' @export
#' 
#' @examples 
#' \dontrun{
#' res <- arg1 %//% 1 # default value if arg1 is not defined
#' }
#' @rdname operator_slashslash
#' @usage a \%//\% b
"%//%" <- function(a, b) {
  if (is.null(a)) b else a
}

#' set an on.exit handler in another function/env
#'
#' on.exit can normally only be set in the current function.
#' Trying to set it using eval() or delayed evaluation does not work
#'
#' See this dicussion:
#' \url{https://stat.ethz.ch/pipermail/r-devel/2013-November/067874.html}
#'
#' @param expr    the expression to be executed by on.exit. the expression
#'                must be evaluable in the where environment
#' @param where   in which function environment to execute/set the on.exit
#'
#' @export
#' @seealso \code{\link{setup_temp_dir}}
add_on_exit <- function(expr, where = parent.frame()) {
  do.call("on.exit", list(substitute(expr), add = TRUE), envir = where)
}


