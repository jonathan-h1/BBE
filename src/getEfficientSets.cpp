#include <Rcpp.h>
#include <math.h>
#include <algorithm>
#include <queue>
using namespace Rcpp;

struct domComp_t {
 domComp_t(std::vector< std::vector<int> >& rks, int nrks) : ranks(rks), nRanks(nrks) {}
 bool operator() (long i,long j) {
   for(int ind = 0; ind < nRanks; ind++){
     if(ranks[i][ind] > ranks[j][ind])
       return true;
     if(ranks[i][ind] < ranks[j][ind])
       return false;
   }
   return (true);
   }
 private:
 std::vector< std::vector<int> >& ranks;
 int nRanks;
};

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
int findUnhandeledPt(std::vector<int> handledPts, int start = 0){
  int nPoints = handledPts.size();
  int unhandeledPt = -1;
  for(int i = start; i < nPoints; i++){
    if(handledPts[i] == -1){
      unhandeledPt = i;
      break;
    }
  }
  return unhandeledPt;
}

std::vector<int> findNeighbors(
    std::vector< std::vector<int> > efficientPoints, std::vector<int> point,
    std::vector<int> handledPts, int start, int nPoints, int nDim)
{
  std::vector<int> neighbors;
  std::vector<int> curPt;
  bool isNeighbor;

  for(int i = start; i < nPoints; i++){
    if(handledPts[i] != -1){
      continue;
    }
    curPt = efficientPoints[i];
    isNeighbor = true;
    for(int d = nDim - 1; d >= 0; d--){
      if(abs(curPt[d] - point[d]) > 1){
        isNeighbor = false;
        break;
      }
    }
    if(isNeighbor){
      neighbors.push_back(i);
    }
    if((curPt[nDim - 1] - point[nDim - 1]) > 1){
      break;
    }
  }
  return neighbors;
}

//' For a vector of indices with efficient points within a grid identify efficient sets.
//'
//' @param efficientPoints The vector containing the indices of all efficient points.
//' @param gridSize The side length of the grid.
//' @param nDim The number of dimensions in the decision space.
//' @export
// [[Rcpp::export]]
List getEfficientSets(
    NumericVector efficientPoints,
    int gridSize,
    int nDim,
    bool domSort = false,
    NumericVector rank = NumericVector::create(),
    int nRank = 0
){
  List efficientSets = List::create();
  List tmpEfficientSets = List::create();
  int nPoints = efficientPoints.length();
  std::vector< std::vector<int> > formEffPoints(nPoints, std::vector<int>(nDim,0));
  std::vector<int> neighbors;
  std::vector<int> handledPts(nPoints, -1);
  std::queue<int> ptQueue;
  int curPtInd;
  int unhandeledPt = 0;
  NumericVector efficientSet;

  long rem;
  for(int i = 0; i < nPoints; i++){
    rem = (long) efficientPoints[i];
    for(int d = nDim - 1; d >= 0; d--){
      formEffPoints[i][d] = (int) floor(rem / pow(gridSize, d));
      rem -= (long) formEffPoints[i][d] * pow(gridSize, d);
    }
  }

  while(unhandeledPt != -1){
    handledPts[unhandeledPt] = unhandeledPt;
    ptQueue.push(unhandeledPt);

    while(!ptQueue.empty()){
      curPtInd = ptQueue.front();
      ptQueue.pop();
      efficientSet.push_back(curPtInd);
      neighbors = findNeighbors(formEffPoints, formEffPoints[curPtInd], handledPts, unhandeledPt, nPoints, nDim);
      for(int neighbor : neighbors){
        ptQueue.push(neighbor);
        handledPts[neighbor] = neighbor;
      }

    }



    efficientSets.push_back(efficientSet);
    efficientSet = NumericVector::create();
    unhandeledPt = findUnhandeledPt(handledPts, unhandeledPt);
  }

  int nEffSets = efficientSets.length();

  if(domSort){
    tmpEfficientSets = efficientSets;
    efficientSets = List::create();
    std::vector<int> sortVec(nEffSets, 0);
    std::vector< std::vector<int> > effSetRanks(nEffSets, std::vector<int>(nRank, 0));
    for(int i = 0; i < nEffSets; i++){
      // Rcout << i << ": ";
      NumericVector effPoints = tmpEfficientSets[i];
      for(auto effPoint : effPoints){
        effSetRanks[i][rank[effPoint] - 1]++;
        // Rcout << rank[effPoint] << ", ";
      }
      // Rcout << "\n";
      sortVec[i] = i;
    }

    struct domComp_t domComp(effSetRanks, nRank);
    std::sort(sortVec.begin(), sortVec.end(), domComp);

    for(int i: sortVec){
      Rcout << i << ", ";
      efficientSets.push_back(tmpEfficientSets[i]);
    }
    Rcout << std::endl;

  }

  for(int i = 0; i < nEffSets; i++){
    // Rcout << i << ": ";
    NumericVector effPoints = efficientSets[i];
    for(int j = 0; j < effPoints.length(); j++){
      effPoints[j] = efficientPoints[effPoints[j]];
      // Rcout << rank[effPoint] << ", ";
    }
  }

  return efficientSets;
}

// create list of points (represented as vektor)
// create vector for handeled points
// fifo queue for first basin
// start with point
// search for unhandeled neighbors (break if last dim greater last dim + 1)
// add to queue
// mark points in handled vektor
// add set to list start with next


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

