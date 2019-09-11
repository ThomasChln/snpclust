
context('utils')

test_sample_impute <- function() {
  num_mat <- matrix(1:16, 4)
  num_mat[2:3, ] <- matrix(rep(NA, 8), 2)
  num_mat <- snpclust:::sample_impute(num_mat)
  apply(num_mat, 2, function(column) {
      expect_true(all(column[2:3] %in% column[c(1,4)]))
    })
  char_mat <- matrix(letters[1:16], 4)
  char_mat[2:3, ] <- matrix(rep(NA, 8), 2)
  char_mat <- snpclust:::sample_impute(char_mat)
  apply(char_mat, 2, function(column) {
      expect_true(all(column[2:3] %in% column[c(1,4)]))
    })
}
test_that('sample_impute', test_sample_impute())

test_merge_dfs <- function() {
  dfs <- lapply(1:3, function(i) iris[sample(1:150, 100), ])
  expected <- merge(merge(dfs[[1]], dfs[[2]]), dfs[[3]])
  expect_identical(expected, merge_dfs(dfs))
}
test_that('merge_dfs', test_merge_dfs())

.df_rbind_all <- function() {
  one <- mtcars[1:4, ]
  two <- mtcars[11:14, ]

  expected <- mtcars[c(1:4, 11:14), ]
  rownames(expected) <- NULL

  expect_identical(df_rbind_all(one, two), expected)
  expect_is(df_rbind_all(one, two, one), "data.frame")
  expect_identical(df_rbind_all(list(one, two)), expected)


  expect_identical(df_rbind_all(one, two, use_row_names = TRUE),
     mtcars[c(1:4, 11:14), ])

  expect_identical(df_rbind_all(one, two, use_row_names = TRUE),
     mtcars[c(1:4, 11:14), ])

  ### factors
  df <- df_rbind_all(iris, iris)
  expect_true(is.factor(df$Species))

  df1 <- data.frame(x = as.factor(letters[1:2]))
  df <- df_rbind_all(df1, df1)
  expect_true(is.factor(df$x))

  ### different levels
  df2 <- data.frame(x = as.factor(letters[3:4]))
  df <- df_rbind_all(df1, df2)
  expect_false(is.factor(df$x))

  ### included levels
  df3 <- data.frame(x = as.factor(letters[2]))
  df <- df_rbind_all(df1, df3)
  expect_false(is.factor(df$x))
}
test_that('df_rbind_all', .df_rbind_all())

