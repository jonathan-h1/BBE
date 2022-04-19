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
grid_size =  300
nDim = 2

design <- moPLOT::generateDesign(fn, points.per.dimension = grid_size)
design$obj.space <- moPLOT::calculateObjectiveValues(design$dec.space, fn)


gradients <- moPLOT::computeGradientFieldGrid(design)
divergence <- moPLOT::computeDivergenceGrid(gradients$multi.objective, design$dims, design$step.sizes)
less <- moPLOT::localEfficientSetSkeleton(design, gradients, divergence, integration = "fast")


nonDomSort <- ecr::doNondominatedSorting(t(design$obj.space[less$sinks, ]))
design$efficientSets <- getEfficientSets(less$sinks, grid_size, nDim,
                                         domSort = TRUE, nonDomSort$ranks, length(unique(nonDomSort$ranks)))

cumPathlength <- moPLOT::computeCumulatedPathLengths(design$dec.space, gradients$multi.objective, less$sinks)


test_that("getBasinLabels returns the correct labels", {
  expect_snapshot_value(getBasinLabelsCPP(design$efficientSets, cumPathlength$last.visited), style = 'json2')
})
