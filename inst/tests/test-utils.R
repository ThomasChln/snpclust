source('helper_setup.R')
context('utils') 

.setup_temp_dir <- function() {
  previous <- getwd()

  ### call it from a function with chdir=TRUE
  dir <- (function() {
      hello <- 'world'
      dir <- setup_temp_dir(TRUE)
      expect_true(file.exists(dir))
      expect_equal(getwd(), dir)

      saveRDS(1, 'toto.rds')

      dir
    })()

  expect_equal(getwd(), previous)
  expect_false(file.exists(dir))

  ### call it from a function with chdir=FALSE
  dir <- (function() {
      hello <- 'world'
      dir <- setup_temp_dir(FALSE)

      expect_true(file.exists(dir))
      expect_equal(getwd(), previous)

      saveRDS(1, file.path(dir, 'toto.rds'))

      dir
    })()

  expect_equal(getwd(), previous)
  expect_false(file.exists(dir))

}
test_that('setup_temp_dir', .setup_temp_dir())

..build_error_message <- function() {
  mkmsg <- qbutils:::.build_error_message

  expect_error(mkmsg('ints=%i', 'coucou'), 'invalid format')
  expect_error(mkmsg('ints=%i', 1:3), 'invalid format')

  msg <- mkmsg('no fmt')
  expect_identical(msg, 'no fmt')

  msg <- mkmsg('ints=%s', 1:3)
  expect_identical(msg, 'ints=1,2,3')

  msg <- mkmsg('chars=%s', LETTERS[1:5])
  expect_identical(msg, "chars=A,B,C,D,E")

  msg <- mkmsg("i = %i s = %s",3,'toto')
  expect_identical(msg, 'i = 3 s = toto')

}
test_that('build_error_message', ..build_error_message())


.die_unless <- function() {
  # bad args
  expect_error(die_unless(1), 'Bad logical argument')

  expect_error(die_unless(TRUE), 'missing')
	expect_error(die_unless(FALSE, 'argh'), 'argh')
	expect_error(die_if(TRUE, ' '))
	expect_identical(die_unless(TRUE,' '), NULL)
	expect_identical(die_if(FALSE,' '), NULL)

  # logical vectors: act as any(!cond)
  die_unless(c(TRUE, TRUE), 'argh')
  expect_error(die_unless(c(TRUE, FALSE), 'argh'), 'argh')

  expect_error(die_if(c(TRUE, FALSE), 'argh'), 'argh')
  expect_error(die_if(c(TRUE, TRUE), 'argh'), 'argh')
  expect_that(die_if(c(FALSE, FALSE), 'argh'), not(throws_error()))

	# messages
	expect_error(die_unless(FALSE, "i = %i s = %s",3,'toto'),
    regexp = "i = 3 s = toto")
	expect_error(die_if(TRUE, "%s, %s, %s",'un', 'deux', 'trois'),
    regexp = "un, deux, trois")

  # test vectors in additional arguments (...) for sprintf
  expect_error(die_unless(FALSE, 'ints=%s', 1:3), 'ints=1,2,3')

  # empty condition (logical)
  expect_error(die_if(NULL == 1, 'Argh'), 'Bad logical argument')
  expect_error(die_unless(NULL == 1, 'Argh'), 'Bad logical argument')

}
test_that('die_unless', .die_unless())

.check_fx_args <- function() {
  ### edge cases
  # empty specs
  f <- function(x) {
    check_fx_args()
    TRUE
  }
  expect_error(f(1), 'empty')

  # bad spec arg names
  f <- function(x) {
    check_fx_args(y = 'i*')
    TRUE
  }
  expect_error(f(1), 'unknown')
  # bad spec format
  f <- function(x) {
    check_fx_args(y = 'Z')
    TRUE
  }
  expect_error(f(1), 'spec')


  ### normal cases
  f <- function(x, y, z) {
    check_fx_args(x = 'i*', y = 'I2', z = 's+')
    TRUE
  }

  expect_true(f(1L, 1:2, 'coucou'))

  expect_error(f(1.2, 1, 'coucou'), regexp = 'x')
  expect_error(f(1, 1:3, 'coucou'), regexp = 'y')
  expect_error(f(1, 1:2, ''), regexp = 'y')

}
test_that('check_fx_args', .check_fx_args())



.parseArgSpec <- function() {
  p <- snpclust:::.parse_arg_spec

  # bad formats
  expect_error(p(NULL))
  expect_error(p(''))
  expect_error(p('i'))
  expect_error(p('n!'))
  expect_error(p('Ce'))
  expect_error(p('l*?'))
  expect_error(p('e+'))

  expect_equal(p('i*'),
    list(type = 'integer', length = '*', na = TRUE, null = TRUE))
  expect_equal(p('i?'),
    list(type = 'integer', length = 0:1, na = TRUE, null = TRUE))
  expect_equal(p('!n+'),
    list(type = 'numeric', length = '+', na = TRUE, null = FALSE))
  expect_equal(p('C5'),
    list(type = 'character', length = 5, na = FALSE, null = TRUE))
  expect_equal(p('A3'),
    list(type = 'any', length = 3, na = FALSE, null = TRUE))
  expect_equal(p('S*'),
    list(type = 'string', length = '*', na = FALSE, null = TRUE))
  expect_equal(p('b+'),
    list(type = 'logical', length = '+', na = TRUE, null = TRUE))
  expect_equal(p('l0'),
    list(type = 'list', length = 0, na = TRUE, null = TRUE))
  expect_equal(p('!D1'),
    list(type = 'data.frame', length = 1, na = FALSE, null = FALSE))

  expect_equal(p('ic(5,0)'),
    list(type = 'integer', length = c(0L,5L), na = TRUE, null = TRUE))
  expect_equal(p('i1:10'),
    list(type = 'integer', length = 1:10, na = TRUE, null = TRUE))

  expect_equal(p('i3+'),
    list(type = 'integer', length = '3+', na = TRUE, null = TRUE))
  expect_equal(p('i33+'),
    list(type = 'integer', length = '33+', na = TRUE, null = TRUE))
  expect_error(p('i3.3+'), 'Bad spec: only allowed')
  expect_error(p('i+3'), 'Bad spec: only allowed')
}
test_that('.parse_arg_spec', .parseArgSpec())



.check_arg <- function() {
  # test bad args
  expect_error(check_arg(1, 1), 'Bad spec')
  expect_error(check_arg('i*', 1.21, 1), 'Bad argname')
  expect_error(check_arg('i*', 1.21, 'x', 1), 'Bad call')

  ### integers
  expect_error(check_arg('i*', 1.21), 'Error checking arg')
  expect_error(check_arg('i1', 'coucou'), 'Error checking arg')
  expect_error(check_arg('i+', integer(), 'Error checking arg'))
  expect_error(check_arg('i5', 1:4), 'Error checking arg')
  expect_error(check_arg('i?', 1:4), 'Error checking arg')
  expect_error(check_arg('!I+', NULL), 'Error checking arg')
  expect_error(check_arg('I+', c(1, NA, 2), 'Error checking arg'))
  expect_error(check_arg('!i*', NULL), 'Error checking arg')

  expect_true(check_arg('i*', integer()))
  expect_true(check_arg('i*', NULL))
  expect_true(check_arg('i?', 1L))
  expect_true(check_arg('I2', c(5L, 0L)))
  expect_true(check_arg('I+', NULL))

  ### numeric
  expect_error(check_arg('n1', 'coucou'), 'Error checking arg')
  expect_error(check_arg('n+', numeric(), 'Error checking arg'))
  expect_error(check_arg('!n1', NULL), 'Error checking arg')
  expect_error(check_arg('n2', 1), 'Error checking arg')
  expect_error(check_arg('N+', c(1, NA, 2), 'Error checking arg'))
  expect_error(check_arg('!n*', NULL), 'Error checking arg')

  expect_true(check_arg('n1', NULL))
  expect_true(check_arg('n*', NULL))
  expect_true(check_arg('N3', 1:3 + pi))
  expect_true(check_arg('n+', c(1, NA, 2)))
  expect_true(check_arg('n?', NA_real_))

  ### character
  expect_error(check_arg('c1', 1), 'Error checking arg')
  expect_error(check_arg('c+', character(), 'Error checking arg'))
  expect_error(check_arg('!c1', NULL), 'Error checking arg')
  expect_error(check_arg('c2', 'coucou'), 'Error checking arg')
  expect_error(check_arg('C+', c('toto', NA), 'Error checking arg'))
  expect_error(check_arg('!c*', NULL), 'Error checking arg')

  expect_true(check_arg('c1', NULL))
  expect_true(check_arg('c*', NULL))
  expect_true(check_arg('C3', LETTERS[1:3]))
  expect_true(check_arg('c+', c('toto', NA)))
  expect_true(check_arg('C?', NULL))

  ### string
  expect_error(check_arg('s1', 1), 'Error checking arg')
  expect_error(check_arg('s+', ''), 'Error checking arg')
  expect_error(check_arg('s+', character(), 'Error checking arg'))
  expect_error(check_arg('!s1', NULL), 'Error checking arg')
  expect_error(check_arg('s2', 'coucou'), 'Error checking arg')
  expect_error(check_arg('S+', c('toto', NA), 'Error checking arg'))
  expect_error(check_arg('s+', c('', "coucou"), 'Error checking arg'))
  expect_error(check_arg('!S*', NULL), 'Error checking arg')

  expect_true(check_arg('s1', NULL))
  expect_true(check_arg('s*', NULL))
  expect_true(check_arg('s3', LETTERS[1:3]))
  expect_true(check_arg('s+', c('toto', NA)))
  expect_true(check_arg('S?', 'toto'))

  ### logical
  expect_error(check_arg('b1', 1), 'Error checking arg')
  expect_error(check_arg('b+', logical(), 'Error checking arg'))
  expect_error(check_arg('!b1', NULL), 'Error checking arg')
  expect_error(check_arg('b2', TRUE), 'Error checking arg')
  expect_error(check_arg('B+', c(FALSE, NA), 'Error checking arg'))
  expect_error(check_arg('!b*', NULL), 'Error checking arg')

  expect_true(check_arg('b1', NULL))
  expect_true(check_arg('b*', NULL))
  expect_true(check_arg('b2', c(TRUE, FALSE)))
  expect_true(check_arg('b+', c(TRUE, NA)))
  expect_true(check_arg('B?', FALSE))

  ### list
  expect_error(check_arg('l1', 1), 'Error checking arg')
  expect_error(check_arg('l+', list(), 'Error checking arg'))
  expect_error(check_arg('!l1', NULL), 'Error checking arg')
  expect_error(check_arg('l2', list(1), 'Error checking arg'))
  expect_error(check_arg('L+', list(1, NA), 'Error checking arg'))
  expect_error(check_arg('!l*', NULL), 'Error checking arg')

  expect_true(check_arg('l*', NULL))
  expect_true(check_arg('l2', list(1,2)))
  expect_true(check_arg('l+', list(1, NA)))
  expect_true(check_arg('l?', list()))

  ### data frame
  expect_error(check_arg('d*', list(1), 'Error checking arg'))
  expect_error(check_arg('d+', data.frame(), 'Error checking arg'))
  expect_error(check_arg('!d1', NULL), 'Error checking arg')
  expect_error(check_arg('d2', data.frame(i = 1), 'Error checking arg'))
  expect_error(check_arg('!d*', NULL), 'Error checking arg')
  df <- mtcars
  df[3,1] <- NA
  expect_error(check_arg('D+', df), 'Error checking arg')

  expect_true(check_arg('d*', NULL))
  expect_true(check_arg('d2', data.frame(i = 1, j = 2)))
  expect_true(check_arg('d+', df))
  expect_true(check_arg('d?', data.frame(i = 1:10)))

  ### any
  expect_error(check_arg('a1', 1:2), 'Error checking arg')
  expect_error(check_arg('a+', character(), 'Error checking arg'))
  expect_error(check_arg('!a1', NULL), 'Error checking arg')
  expect_error(check_arg('a2', 'coucou'), 'Error checking arg')
  expect_error(check_arg('A+', c('toto', NA), 'Error checking arg'))
  expect_error(check_arg('!a*', NULL), 'Error checking arg')

  expect_true(check_arg('a*', NULL))
  expect_true(check_arg('A3', 1:3))
  expect_true(check_arg('A+', mtcars))
  expect_true(check_arg('a+', c('toto', NA)))
  expect_true(check_arg('A?', NULL))

  #test "N+" argument
  expect_true(check_arg('I3+', 1L:4L))
  expect_true(check_arg('I3+', 1L:3L))
  expect_error(check_arg('I33+', 1L:3L), 'must be > 32')
  expect_error(check_arg('I3+', 1L), 'must be > 2')
  expect_error(check_arg('I+3', 1L), 'Bad spec: only allowed')
  expect_true(check_arg('d1+', data.frame(i = 1, j = 2)))
  expect_true(check_arg('l1+', list(1,2)))
  expect_true(check_arg('b1+', c(TRUE, FALSE)))
  expect_true(check_arg('n2+', c(1, NA, 2)))
  expect_true(check_arg('c2+', c('toto', NA)))
  expect_true(check_arg('s1+', c('toto', NA)))
}
test_that('check_arg', .check_arg ())


test_df_columns_exist <- function (){
  expect_error(df_columns_exist(1), 'cols')
  expect_error(df_columns_exist('dummy', 1), 'df')

  #one colname
  df <- data.frame(toto = 2, titi = 3)
  expect_true(df_columns_exist(df, 'toto', silent = TRUE))
  expect_false(df_columns_exist(df, 'dummy', silent = TRUE))

  #more than one colname
  expect_false(df_columns_exist(df, c('toto', 'dummy'), silent = TRUE))
  expect_true(df_columns_exist(df, c('toto', 'titi'), silent = TRUE))
  expect_true(df_columns_exist(df, c('titi', 'toto'), silent = TRUE))  #order doesn't natter

  #test warning for one colname not in df
  expect_warning(df_columns_exist(df, c('toto', 'dummy')), 'dummy')
  #test warning for two colnames not in df
  expect_warning(df_columns_exist(df, c('tutu', 'dummy')), 'dummy, tutu')
}
test_that('df_columns_exist',test_df_columns_exist())


.catch_warnings <- function() {

  # do not catch errors
  expect_error(catch_warnings(stop('argh')))

  ll <- catch_warnings(1+1)
  expect_equal(ll[[1]], 2)
  expect_equal(ll[[2]], list())

  ll <- catch_warnings({warning('Argh'); iris  })
  expect_equal(ll[[1]], iris)
  expect_equal(ll[[2]][[1]]$message, 'Argh')

  ll <- catch_warnings({warning('w1'); warning('w2');
      (function() { warning('w3')})(); pi  })
  expect_equal(ll[[1]], pi)
  ws <- ll$warnings
  expect_true(length(ws)==3 && ws[[1]]$message=='w1' &&
      ws[[2]]$message=='w2' && ws[[3]]$message=='w3')

}
test_that("catch_warnings", .catch_warnings())


