#include <Rcpp.h>
#include <math.h>
#include <algorithm>
#include <queue>
#include <limits>
using namespace Rcpp;

std::vector<long> getNeighbors(long point, int gridSize, int nDim){

  std::vector<long> neighbors;
  neighbors.reserve((int) pow(3, nDim));

  if(nDim == 1){
    for(long i = -1; i < 2; i++){
      if(point + i < gridSize  && point + i >= 0){
        neighbors.push_back(point + i);
      }
    }
    return neighbors;
  }

  long dimStep = (long) pow(gridSize, nDim - 1);
  int dimInd = (int) floor(point / dimStep);
  long rem = point - dimInd * dimStep;
  std::vector<long> prevNeighbors = getNeighbors(rem, gridSize, nDim - 1);
  for(long i = -1; i < 2; i++){
    if(dimInd + i <= gridSize  && dimInd + i >= 0){
      for(long prevPoint: prevNeighbors){
        neighbors.push_back((dimInd + i) * dimStep + prevPoint);
      }
    }
  }
  return neighbors;


}

void handleQueue(
    std::queue<long>& ptQueue, NumericVector& basinLabels, NumericVector& gradients,
    int label, int gridSize, int nDim){

  long curPoint;
  std::vector<long> neighbors;

  while(!ptQueue.empty()){
    curPoint = ptQueue.front();
    ptQueue.pop();
    //std::cout<< "Current Point is: " << curPoint<< std::endl << "found neighbors: "<<std::endl;
    neighbors = getNeighbors(curPoint, gridSize, nDim);
    for(long neighbor: neighbors){
      //std::cout<<'\t'<< neighbor;
      if(basinLabels[neighbor] != -1 || gradients[neighbor] < gradients[curPoint]){
        //std::cout<<std::endl;
        continue;
      }
      //std::cout<<" - labeled"<< std::endl;
      basinLabels[neighbor] = label;
      ptQueue.push(neighbor);
    }
    //std::cout<<std::endl;
  }
}

/*struct gradComp_t {
 gradComp_t(std::vector<int>& grads) : gradients(grads) {}
 bool operator() (long i,long j) { return (gradients[i] < gradients[j]);}
 private:
 std::vector<int>& gradients;
};*/

std::vector<long> getUnhandeled(NumericVector& basinLabels, NumericVector& gradients, long nPoints, int gridSize, int nDim){
  std::vector<long> neighbors;
  std::vector<long> unhandeledPt;
  long ptwSmallestNeigh;
  double smallestNeighGrad = std::numeric_limits<double>::max();
  long ptSmallestNeigh;

  for(long i = 0; i < nPoints; i++){
    if(basinLabels[i] == -1){
      neighbors = getNeighbors(i, gridSize, nDim);
      for(long neighbor: neighbors){
        if(basinLabels[neighbor] != -1 && gradients[neighbor] < smallestNeighGrad){
          smallestNeighGrad = gradients[neighbor];
          ptwSmallestNeigh = i;
          ptSmallestNeigh = neighbor;
        }
      }
    }
  }

  if(smallestNeighGrad == std::numeric_limits<double>::max()){
    return unhandeledPt = {-1, -1};
  }


  unhandeledPt = {ptwSmallestNeigh, ptSmallestNeigh};
  return unhandeledPt;

}

//' Identify for every point in the rasterised decision space the corresponding basin it belongs to.
//'
//' @param efficientSets A list containing for every efficient set the indices of all efficient points.
//' @param gradients A vector containing the accumulated gradients for every point in the grid.
//' @param gridSize The side length of the grid.
//' @param nDim The number of dimensions in the decision space.
//' @export
// [[Rcpp::export]]
NumericVector getBasinLabels(List efficientSets, NumericVector gradients, int gridSize, int nDim){

  long nPoints = (long) pow(gridSize, nDim);
  NumericVector basinLabels (nPoints, -1);
  int nSets = efficientSets.size();
  NumericVector curSet;
  std::vector<int> neighbors;
  std::vector<int> handledPts(nPoints, -1);
  std::queue<long> ptQueue;
  std::vector<long> unhandeledPt;

  for(int set = 1; set <= nSets; set++){
    curSet = efficientSets[set - 1];
    for(auto effPoint : curSet){
      ptQueue.push((long) effPoint);
      basinLabels[(long) effPoint] = set;
    }
    handleQueue(ptQueue, basinLabels, gradients, set, gridSize, nDim);
  }
  while(true){
    unhandeledPt = getUnhandeled(basinLabels, gradients, nPoints, gridSize, nDim);
    if(unhandeledPt[0] == (long) -1){
      break;
    }
    ptQueue.push(unhandeledPt[0]);
    basinLabels[unhandeledPt[0]] = basinLabels[unhandeledPt[1]];
    handleQueue(ptQueue, basinLabels, gradients, basinLabels[unhandeledPt[1]], gridSize, nDim);
  }
  return basinLabels;
}

