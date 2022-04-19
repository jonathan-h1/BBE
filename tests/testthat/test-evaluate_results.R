library(bbe)
library(tidyverse)

fn = smoof::makeDTLZ1Function(2,2)
tb <- tibble(fun_calls = c(rep(10, 11), rep(20, 11)), x1 = c(0:10 * 0.1, 0:10 * 0.1), x2 = c(rep(0.4, 11), rep(0.5, 11)), y1 = rep(0.0, 22), y2 = rep(0.0, 22))
tb[,4:5] <- t(apply(tb, 1, function(x){
  fn(x[2:3])
}))
design <- evaluate_results(tb, fn, ref.point = smoof::getRefPoint(fn))

mean_vec <- c(ecr::computeHV(t(tb[1:11,4:5]), c(11,11)) / 3, ecr::computeHV(t(tb[12:22,4:5]), c(11,11)) / 3)

test_that(
  'The values are caculated correctly when using HV.',
{
  expect_equal(design$basin_separated_eval$value_basin1, c(0.0, ecr::computeHV(t(tb[12:22,4:5]), c(11,11))))
  expect_equal(design$basin_separated_eval$value_basin2, c(ecr::computeHV(t(tb[1:11,4:5]), c(11,11)), 0.0))
  expect_equal(design$basin_separated_eval$value_basin3, c(0.0, 0.0))
  expect_equal(design$basin_separated_eval$mean_value, mean_vec)
  expect_equal(design$basin_separated_eval$auc_hv_mean, c(mean_vec[1] * 10 / 2,
                                                          ((mean_vec[1] + mean_vec[2] - mean_vec[1]) / 2 + mean_vec[1]) * 10))
  expect_equal(design$basin_separated_eval$auc_hv1, c(0.0, ecr::computeHV(t(tb[12:22,4:5]), c(11,11)) * 10 / 2))
}
)
