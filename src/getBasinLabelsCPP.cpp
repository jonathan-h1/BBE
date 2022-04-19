#include <Rcpp.h>
#include <math.h>
#include <algorithm>
#include <queue>
#include <limits>
using namespace Rcpp;

// Identify for every point in the rasterized decision space the corresponding basin it belongs to.
//
// @param efficientSets A list containing for every efficient set the indices of all efficient points.
// @param gradients A vector containing the accumulated gradients for every point in the grid.
// @param gridSize The side length of the grid.
// @param nDim The number of dimensions in the decision space.
// [[Rcpp::export]]
NumericVector getBasinLabelsCPP(List efficientSets, NumericVector lastVisited){

  Rcout << "Computing basin labels ... \n";

  long nPoints = (long) lastVisited.size();
  NumericVector basinLabels (nPoints, -1);
  int nSets = efficientSets.size();
  NumericVector curSet;

  for(int set = 1; set <= nSets; set++){
    curSet = efficientSets(set - 1);
    for(auto effPoint : curSet){
      basinLabels((long) effPoint - 1) = set;
    }
  }

  for(long point = 0; point < nPoints; point++){
    if(basinLabels(point) == -1){
      if(lastVisited(point) == -1){
        continue;
      } else {
        basinLabels(point) = basinLabels(lastVisited(point) - 1);
      }
    }
  }

  return basinLabels;
}


