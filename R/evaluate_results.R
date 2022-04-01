#' Evaluate the performance of an algorithm in a basin separated manner.
#' Internally the basins are determined with the efficient points
#' and the gradients calculated by ['moPlot'] for the rasterized decision space.
#' See the paper for more details. In the following the number of dimensions in
#' the decision space will be denoted as dec.nDim and the number of dimensions in
#' the objective space will be denoted as obj.nDim.
#'
#' @param results \[\code{tibble}\] \cr
#'  A tibble with the results of an algorithm run. It should be
#'  organized as follows: The first column should contain the function calls needed
#'  to retrieve the solutions captured in the remaining columns. The second column
#'  until the 1 + dec.nDim column should contain the coordinates of the solutions
#'  in the decision space. Finally, the remaining columns should contain the
#'  coordinates in the objective space. Note that the point will be grouped by
#'  the function calls needed to retrieve them. Thus, make sure that points that
#'  should be evaluated together have the same value in the first column.
#' @param fn \[\code{function}\] \cr
#'  The multi-objective function under consideration.
#'  It should be a ['smoof']-function. The number of objectives and the upper and
#'  lower bounds are inferred from the function.
#' @param ... \[\code{any}\] \cr
#'  Further arguments that should be passed to ['eval_fn'].
#'  In the default case \code{ref.point} should be passed here
#'  for the calculation of the hypervolume.
#' @param eval_fn \[\code{function}\] \cr
#'  The function that is used to evaluate the solutions in a basin.
#'  It should accept the points in a column wise dataframe. The default is
#'  ['computeHV'].
#' @param grid_size \[\code{integer(1)}\] \cr
#'  The granuality of the raster per dimension. The default is 300.
#' @param basins \[\code{integer}\] \cr
#'  A vector of integers identifying the basins that should be
#'  considered during the evaluation. The default is to consider the first three
#'  basins.
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
#' @param dec_space_labels \[\code{integer} | \code{NULL}\] \cr
#'  If a custom labeling of the grid in the decision space is wanted
#'  the labels can be supplied here. The vector should contain the wanted label
#'  at the position of the index of a cell in the grid. If \code{NULL} the labels
#'  are computed from the points in the efficient sets.
#' @return \[\code{list}\] \cr
#'  A \code{design} list from ['moPLOT']. Additionally attached are
#'  the efficient sets (\code{efficientSets}), the labels for the decision space
#'  (\code{decSpaceLabels}) and a tibble (\code{basin_separated_eval})
#'  with the basin separated results.
#' @examples
#' # NOT RUN {
#' fn <- smoof::makeDTLZ1function(2,2)
#' # dummy tibble
#' tb <- tibble::tibble(fc = 5, x1 = 0.5, x2 = 0.5, y1 = 0.25, y2 = 0.25)
#' evalutate_results(tb, fn, ref.point = smoof::getRefPoint(fn))
#' # }
#' @export
evalutate_results = function(results, fn, ...,
                             eval_fn = ecr::computeHV,
                             grid_size = 300L,
                             basins = 1:3,
                             join_fronts = FALSE,
                             keep_points = FALSE,
                             efficient_sets = NULL,
                             dec_space_labels = NULL,
                             design = NULL){

  assert_data_frame(results, types = c('numeric'), min.cols = 5)
  assert_function(fn)
  assert_function(eval_fn)
  assert_class(fn, c('smoof_function')) #, 'smoof_multi_objective_function' # wrappers and loggers remove this class...
  if(identical(ecr::computeHV, eval_fn)){
    arguments <- list(...)
    assert('ref.point' %in% names(arguments), .var.name = 'ref.point')
  }
  assert_atomic_vector(basins, any.missing = FALSE, min.len = 1, unique = TRUE)

  nDim <- as.integer(smoof::getNumberOfParameters(fn))
  nDimObj <- as.integer(smoof::getNumberOfObjectives(fn))

  if(!is.null(dec_space_labels)){
    assert_atomic_vector(dec_space_labels, len = grid_size ** nDim, any.missing = FALSE)
  }

  if(!is.null(efficient_sets)){
    assert_list(efficient_sets, min.len = 1, any.missing = FALSE)
  }

  if(is.null(design)){
    design <- moPLOT::generateDesign(fn, points.per.dimension = grid_size)
    design$obj.space <- moPLOT::calculateObjectiveValues(design$dec.space, fn)
  }

  if(is.null(dec_space_labels)){
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
  } else {
    design$decSpaceLabels <- dec_space_labels
  }

  if(keep_points){
    list_it <- lapply(unique(results[[1]]), function(x) {
      return(results[results[[1]] <= x, ])
    })
  } else  {
    list_it <- split(results, factor(results[[1]]))
  }

  boundaries <- c(rbind(getLowerBoxConstraints(fn), getUpperBoxConstraints(fn)))
  cat('Evaluating per basin ...\n')
  res_per_basin <- lapply(list_it, function(df_part) {
    # points <- df_part[, 2:1 + nDim]
    points <- select(df_part, paste0("x",c(1:nDim))) # Enforce right selection of columns
    points$labels <- filterByBasin(points, design$decSpaceLabels, boundaries, grid_size, nDim)

    perf_vals <- sapply(basins, function(x){
      basin_points <- df_part[points$labels == x, ]
      perf_val = 0.0
      if (nrow(basin_points) != 0) {
        perf_val <- eval_fn(t(basin_points[, (2 + nDim):(1 + nDim + nDimObj)]), ...)
      }
      return(perf_val)
    })

    names(perf_vals) <- paste0('value_basin', basins)

    mean_val = mean(perf_vals)

    return(data.frame(
      fun_calls = df_part$fun_calls[nrow(df_part)],
      t(perf_vals),
      mean_value = mean_val
    ))
  })

  design$basin_separated_eval <- tibble::as_tibble(do.call("rbind", res_per_basin))
  auc <- design$basin_separated_eval %>%
    mutate(max_hv_mean = pmax(mean_value, c(0, mean_value)[1:nrow(design$basin_separated_eval)]),
           min_hv_mean = pmin(mean_value, c(0, mean_value)[1:nrow(design$basin_separated_eval)]),
           max_hv1 = pmax(value_basin1, c(0, value_basin1)[1:nrow(design$basin_separated_eval)]),
           min_hv1 = pmin(value_basin1, c(0, value_basin1)[1:nrow(design$basin_separated_eval)]),
           fn_calls_diff = fun_calls - c(0, fun_calls)[1:nrow(design$basin_separated_eval)]) %>%
    transmute(auc_hv_mean = cumsum(fn_calls_diff * (min_hv_mean + (max_hv_mean - min_hv_mean) / 2)),
              auc_hv1 = cumsum(fn_calls_diff * (min_hv1 + (max_hv1 - min_hv1) / 2)))

  # fn_calls = c(1, 3, 4, 5, 7)
  # hv_mean = c(2,4,5,4,2)
  # hv1 = c(5, 2, 2, 2, 4)
  # design = list(basin_separated_eval = tibble(fun_calls = fn_calls, mean_value = hv_mean, value_basin1 = hv1))


  design$basin_separated_eval <- cbind(design$basin_separated_eval, auc)

  return(design)
}
