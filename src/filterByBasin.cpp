#include <Rcpp.h>
#include <math.h>
#include <iostream>
using namespace Rcpp;




//' Group a set of points by the label of the region in the rasterized decision space htey lie in.
//'
//' @param solutions A data frame with one column per dimension.
//' @param basinLabels A vector containing the labels of all regions as retrieved from getBasinLabels.
//' @param boundaries A vector containing the minimal and maximal possible values per dimension.
//' @param nBasins The number of unique labels in [`basinLabels`].
//' @param gridSize The side length of the grid.
//' @param nDim The number of dimensions in the decision space.
// [[Rcpp::export]]
NumericVector filterByBasin(DataFrame solutions, NumericVector basinLabels, NumericVector boudaries,
                   int nBasins, int gridSize, int nDim){

  int nSolutions = solutions.nrow();
  NumericVector solutionLabels(nSolutions);
  // for(int i = 0; i < nBasins; i++){
  //   filteredSolutions.push_back(NumericVector::create());
  // }
  NumericVector dimVals;
  int label;
  std::vector<double> minAndStep;
  minAndStep.reserve(nDim * 2);
  for (int dim = 0; dim < nDim; ++dim)
  {
    minAndStep.push_back((double) boudaries(dim * 2));
    minAndStep.push_back((double) ((boudaries(dim * 2 + 1) - boudaries(dim * 2)) / gridSize));
  }


  for (int i = 0; i < nSolutions; ++i)
  {
    long index = 0;
    for (int dim = 0; dim < nDim; ++dim)
    {
      dimVals = solutions[dim];
      int dInd = (int) floor((dimVals(i) - minAndStep.at(dim * 2)) / minAndStep.at(dim * 2 + 1));
      if(dInd == gridSize){
        dInd--;
      }
      index += (long) dInd * pow(gridSize, dim);
    }
    // std::cout<<"Index: "<< index<<std::endl;

    label = basinLabels(index);
    // std::cout<<"Label: "<<label<<std::endl;
    solutionLabels(i) = label;
  }
  return solutionLabels;
}


