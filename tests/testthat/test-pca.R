context('pca functions')
library(magrittr)

pca = iris %>% cbind(id = rownames(.), .) %>%
  get_pca('id', obs = 1:100 %>% as.character, var = names(iris)[1:3])

test_pca_fortify = function() {
  expect_is(snpclust:::pca_fortify(pca), 'data.frame')
}
test_that('pca_fortify', test_pca_fortify())

test_pca_vars_plot = function() {
  gg = ggplot_pca(pca, obs_sup = TRUE, vars = TRUE, vars_sup = TRUE)
  expect_is(ggplot2::ggplotGrob(gg), 'gtable')
}
test_that('pca_vars_plot', test_pca_vars_plot())
