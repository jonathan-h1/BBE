library(bbe)
library(moPLOT)
library(tidyverse)
library(smoof)

#### Aspar Function ####
f1_1 = function(x) (x[1]**4 - 2*x[1]**2 + x[2]**2 + 1)
f2_1 = function(x) ((x[1] + 0.5)**2 + (x[2]-2)**2)
f3_1 = function(x) ((x[1] + 0.25) ** 4 + 3 * (x[2] - 1) ** 2)
f_2d2d = function(x) c(f1_1(x), f2_1(x))
f_2d3d = function(x) c(f1_1(x), f2_1(x), f3_1(x))

f1_2 = function(x) (x[1]**4 - 2*x[1]**2 + x[2]**2 + 1 + x[3] ** 2)
f2_2 = function(x) ((x[1] + 0.5)**2 + (x[2]-2)**2 + (x[3] - 1) ** 4)
f3_2 = function(x) ((x[1] + 0.25) ** 4 + 3 * (x[2] - 1) ** 2 + (x[3] - 1) ** 2)
f_3d2d = function(x) c(f1_2(x), f2_2(x))
f_3d3d = function(x) c(f1_2(x), f2_2(x), f3_2(x))

makeAsparFunction <- function(dimensions = 2, n.objectives = 2) {
  if (dimensions == 2 && n.objectives == 2) {
    smoof::makeMultiObjectiveFunction(name = "2D->2D Test Function", id = "test_2d2d", description = "", fn = f_2d2d,
                                      par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = c(-2,-1), upper = c(2,3)))
  } else if (dimensions == 2 && n.objectives == 3) {
    smoof::makeMultiObjectiveFunction(name = "2D->3D Test Function", id = "test_2d3d", description = "", fn = f_2d3d,
                                      par.set = ParamHelpers::makeNumericParamSet(len = 2, lower = c(-2,-1), upper = c(2,3)))
  } else if (dimensions == 3 && n.objectives == 2) {
    smoof::makeMultiObjectiveFunction(name = "3D->2D Test Function", id = "test_3d2d", description = "", fn = f_3d2d,
                                      par.set = ParamHelpers::makeNumericParamSet(len = 3, lower = c(-2,-1,-2), upper = c(2,3,2)))
  } else if (dimensions == 3 && n.objectives == 3) {
    smoof::makeMultiObjectiveFunction(name = "3D->3D Test Function", id = "test_3d3d", description = "", fn = f_3d3d,
                                      par.set = ParamHelpers::makeNumericParamSet(len = 3, lower = c(-2,-1,-2), upper = c(2,3,2)))
  }
}
#### end Aspar ####

fn <- makeAsparFunction()

design <- moPLOT::generateDesign(fn, points.per.dimension = 300)
design$obj.space <- calculateObjectiveValues(design$dec.space, fn, parallelize = F)

gradients <- computeGradientFieldGrid(design)
divergence <- computeDivergenceGrid(gradients$multi.objective, design$dims, design$step.sizes)
less <<- localEfficientSetSkeleton(design, gradients, divergence, integration="fast")


test_that("efficient sets are identified with efficient points", {
  expect_equal(getEfficientSets(c(1,2,3,11,13,15, 16,27,38,49,59,69), 10L, 2L,
                                domSort = TRUE,
                                rank = c(1,1,1,1,2,1,1,2,2,2,2,1), nRank = 2),
               list(c(1,2,11,3,13), c(15, 16,27,38,49,59,69)))
  expect_snapshot_value(getEfficientSets(less$sinks, 300L, 2L), style = 'json2')
})
