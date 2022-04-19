#' @description Compute the basin labels for the points in the grid structure
#' for a given function as well as the efficient sets. To be further used by
#' ['evaluate_results'].
#'
#'
#' @param fn \[\code{function}\] \cr
#'  The multi-objective function under consideration.
#'  It should be a ['smoof']-function. The number of objectives and the upper and
#'  lower bounds are inferred from the function.
#' @param grid_size \[\code{integer(1)}\] \cr
#'  The granuality of the raster per dimension. The default is 300.
#' @param join_fronts \[\code{logical(1)}\] \cr
#'  This should be \code{TRUE} if the efficient sets should be joined when
#'  they are part of the same domination front. Default is \code{FALSE}.
#' @param efficient_sets \[\code{list} | \code{NULL}\] \cr
#'  If the efficient points returned by ['moPLOT'] are not accurate or enought and
#'  more precice locations are available they or a custom merging of sets is wanted
#'  this can be supplyed here. Expected is a list of vectors containing the indices
#'  of the gridcells that should be regarded an efficient point.
#'  An element of the List is treated as one efficient set.
#'  If \code{NULL} the efficient points returned by ['moPLOT'] are merged
#'  to efficient sets and those are then ordered by the number of points
#'  in a domination layer. Default is \code{NULL}.
#' @return \[\code{list}\] \cr
#'  A \code{design} list from ['moPLOT']. Additionally attached are
#'  the efficient sets (\code{efficientSets}) and the labels for the decision space
#'  (\code{decSpaceLabels}).
#' @export
get_decision_space_labels = function(fn,
                             grid_size = 300L,
                             join_fronts = FALSE,
                             efficient_sets = NULL){

  checkmate::assert_function(fn)
  checkmate::assert_class(fn, c('smoof_function', 'smoof_multi_objective_function'))

  nDim <- as.integer(smoof::getNumberOfParameters(fn))
  nDimObj <- as.integer(smoof::getNumberOfObjectives(fn))

  if(!is.null(efficient_sets)){
    checkmate::assert_list(efficient_sets, min.len = 1, any.missing = FALSE)
  }

  design <- moPLOT::generateDesign(fn, points.per.dimension = grid_size)
  design$obj.space <- moPLOT::calculateObjectiveValues(design$dec.space, fn)


  gradients <- moPLOT::computeGradientFieldGrid(design)
  divergence <- moPLOT::computeDivergenceGrid(gradients$multi.objective, design$dims, design$step.sizes)
  less <- moPLOT::localEfficientSetSkeleton(design, gradients, divergence, integration = "fast")

  if(is.null(efficient_sets)){
    # print('Computing efficient sets')
    nonDomSort <- ecr::doNondominatedSorting(t(design$obj.space[less$sinks, ]))
    design$efficientSets <- getEfficientSets(less$sinks, grid_size, nDim,
                                             domSort = TRUE, nonDomSort$ranks, length(unique(nonDomSort$ranks)),
                                             joinFronts = join_fronts)
  } else {
    design$efficientSets <- efficient_sets
  }
  cumPathlength <- moPLOT::computeCumulatedPathLengths(design$dec.space, gradients$multi.objective, less$sinks)
  design$decSpaceLabels <- getBasinLabelsCPP(design$efficientSets, cumPathlength$last.visited)

  return(design)
}
