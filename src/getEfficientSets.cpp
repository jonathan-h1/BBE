#include <Rcpp.h>
#include <math.h>
#include <algorithm>
#include <queue>
using namespace Rcpp;

struct domComp_t {
 domComp_t(std::vector< std::vector<int> >& rks, int nrks) : ranks(rks), nRanks(nrks) {}
 bool operator() (long i,long j) {
   for(int ind = 0; ind < nRanks; ind++){
     if(ranks.at(i).at(ind) > ranks.at(j).at(ind))
       return true;
     if(ranks.at(i).at(ind) < ranks.at(j).at(ind))
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
    if(handledPts.at(i) == -1){
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
    if(handledPts.at(i) != -1){
      continue;
    }
    curPt = efficientPoints.at(i);
    isNeighbor = true;
    for(int d = nDim - 1; d >= 0; d--){
      if(abs(curPt.at(d) - point.at(d)) > 1){
        isNeighbor = false;
        break;
      }
    }
    if(isNeighbor){
      neighbors.push_back(i);
    }
    if((curPt.at(nDim - 1) - point.at(nDim - 1)) > 1){
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
// [[Rcpp::export]]
List getEfficientSets(
    NumericVector efficientPoints,
    int gridSize,
    int nDim,
    bool domSort = false,
    NumericVector rank = NumericVector::create(),
    int nRank = 0,
    bool joinFronts = false
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

  Rcout << "Computing efficient sets ... \n";
  // Rcout << "nPoints = " << nPoints << " nIndividualRanks = " << rank.size() << std::endl;

  long rem;
  for(int i = 0; i < nPoints; i++){
    rem = (long) efficientPoints(i) - 1;
    for(int d = nDim - 1; d >= 0; d--){
      formEffPoints.at(i).at(d) = (int) floor(rem / pow(gridSize, d));
      rem -= (long) formEffPoints.at(i).at(d) * pow(gridSize, d);
    }
  }

  while(unhandeledPt != -1){
    handledPts.at(unhandeledPt) = unhandeledPt;
    ptQueue.push(unhandeledPt);

    while(!ptQueue.empty()){
      curPtInd = ptQueue.front();
      ptQueue.pop();
      efficientSet.push_back(curPtInd);
      neighbors = findNeighbors(formEffPoints, formEffPoints.at(curPtInd), handledPts, unhandeledPt, nPoints, nDim);
      for(int neighbor : neighbors){
        ptQueue.push(neighbor);
        handledPts.at(neighbor) = neighbor;
      }

    }



    efficientSets.push_back(efficientSet);
    efficientSet = NumericVector::create();
    unhandeledPt = findUnhandeledPt(handledPts, unhandeledPt);
  }

  int nEffSets = efficientSets.length();
  std::vector<int> sortVec(nEffSets, 0);
  std::vector< std::vector<int> > effSetRanks(nEffSets, std::vector<int>(nRank, 0));

  if(domSort){
    // Rcout << "Dom Sort: neffSets = " << nEffSets << std::endl;
    tmpEfficientSets = efficientSets;
    efficientSets = List::create();
    for(int i = 0; i < nEffSets; i++){
      NumericVector effPoints = tmpEfficientSets(i);
      for(auto effPoint : effPoints){
        effSetRanks.at(i).at(rank(effPoint) - 1)++;
      }
      sortVec.at(i) = i;
    }

    struct domComp_t domComp(effSetRanks, nRank);
    std::sort(sortVec.begin(), sortVec.end(), domComp);

    if(!joinFronts){
      for(int i: sortVec){
        efficientSets.push_back(tmpEfficientSets(i));
      }
    }

  }

  if(joinFronts && domSort){
    // Rcout << "Join Fronts" << std::endl;
    efficientSets = List::create();
    int currentLayer = 0;
    int start = 0;
    std::vector<int> front;
    front.reserve(nRank);
    std::vector<int>::iterator it;

    for(int effSetInd = 0; effSetInd < nEffSets; effSetInd++){
      if(effSetRanks.at(sortVec.at(effSetInd)).at(currentLayer) > 0){
        front.push_back(sortVec.at(effSetInd));
      }
      else if(front.empty()){
        currentLayer++;
        effSetInd--;
      } else {
        currentLayer++;
        effSetInd--;
        std::vector<int> joinedVec;
        joinedVec.reserve((int) ((nPoints / nEffSets) * 2 * front.size()));
        for(int i : front){
          efficientSet = tmpEfficientSets(i);
          joinedVec.insert(joinedVec.begin(), efficientSet.begin(), efficientSet.end());
        }
        front.clear();
        efficientSet.assign(joinedVec.begin(), joinedVec.end());
        efficientSets.push_back(efficientSet);
      }
    }
  }

  for(int i = 0; i < efficientSets.size(); i++){
    NumericVector effPoints = efficientSets(i);
    for(int j = 0; j < effPoints.length(); j++){
      effPoints(j) = efficientPoints(effPoints(j));
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

