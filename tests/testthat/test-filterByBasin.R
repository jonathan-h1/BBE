library(bbe)

data <- data.frame(x = c(1.5, 5.5, 7.9, 10.0), y=c(15, 9, 4, 20))
bounds <- c(0.0, 10.0, 0.0, 20.0)
labels <- c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2)


test_that("the right labels are found for a set of solutions", {
  expect_equal(filterByBasin(data, labels, bounds, 2,4,2)[[3]], c(2,1,1,2))
})
