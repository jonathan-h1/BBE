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
    if(dimInd + i < gridSize  && dimInd + i >= 0){
      for(long prevPoint: prevNeighbors){
        neighbors.push_back((dimInd + i) * dimStep + prevPoint);
      }
    }
  }
  //~pervNeighbors;
  return neighbors;


}

void handleQueue(
    std::queue<long>& ptQueue, NumericVector& basinLabels, NumericVector& gradients,
    std::vector<int>& dists, int label, int gridSize, int nDim,
    bool printLog = false){

  long curPoint;
  std::vector<long> neighbors;
  bool boundary = false;
  int dinInt;
  long rem;

  while(!ptQueue.empty()){
    curPoint = ptQueue.front();
    ptQueue.pop();
    //rem = curPoint;
    /*for(int i = nDim - 1; i >=0; i--){
     dinInt = floor(rem / pow(gridSize, i));
     if(dinInt == 0 || dinInt == (gridSize - 1)){
     boundary = true;
     }
     rem -= dinInt * pow(gridSize, i);
    }
     if(boundary){
     continue;
     }*/
    if(printLog){
      Rcout<< "Pt: " << curPoint<< " (" << label << ") " << "-- {";
    }
    neighbors = getNeighbors(curPoint, gridSize, nDim);
    bool start = true;
    for(long neighbor: neighbors){
      if(printLog && !(start)){
        Rcout<< ", " << neighbor <<"(";
      }
      if(printLog && start){
        Rcout<< neighbor <<"(";
        start = false;
      }
      if(gradients(neighbor) < gradients(curPoint)){
        if(printLog){
          Rcout<< "<)";
        }
        continue;
      }
      // std::cout<< "dist point: " << dists[curPoint] << " + 1 " << (((int) dists[curPoint]) + 1) << " dist neigh " << dists[neighbor] <<std::endl;
      if(basinLabels(neighbor) != -1 && ((((int) dists.at(curPoint)) + 1) >= dists.at(curPoint))){
        if(printLog){
          Rcout<<"-)";
        }
        continue;
      }
      if(printLog){
        Rcout<<"+)";
      }
      basinLabels(neighbor) = label;
      dists.at(neighbor) = ((int) dists.at(curPoint)) + 1;
      ptQueue.push(neighbor);
    }
    if(printLog){
      Rcout<< "}" << std::endl;
    }
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
    if(basinLabels(i) == -1){
      neighbors = getNeighbors(i, gridSize, nDim);
      for(long neighbor: neighbors){
        if(basinLabels(neighbor) != -1 && gradients(neighbor) < smallestNeighGrad){
          smallestNeighGrad = gradients(neighbor);
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
  std::vector<int> dists (nPoints, 0);
  //std::cout<< "labeling sets" << std::endl;
  for(int set = 1; set <= nSets; set++){
    //std::cout<< "Set: " << set << std::endl;
    curSet = efficientSets(set - 1);
    for(auto effPoint : curSet){
      ptQueue.push((long) effPoint - 1);
      basinLabels((long) effPoint - 1) = set;
      dists.at((long) effPoint - 1) = 0;
    }
    handleQueue(ptQueue, basinLabels, gradients, dists, set, gridSize, nDim);
  }
  // for( int i = 0; i < dists.size(); i++){
  //   if(i % gridSize == 0) Rcout << std::endl;
  //   Rcout << dists[i] << "";
  // }
  // Rcout << std::endl;
  //std::cout<< "before while" << std::endl;
  while(true){
    unhandeledPt = getUnhandeled(basinLabels, gradients, nPoints, gridSize, nDim);
    if(unhandeledPt.at(0) == (long) -1){
      break;
    }
    ptQueue.push(unhandeledPt[0]);
    basinLabels(unhandeledPt[0]) = basinLabels(unhandeledPt[1]);
    dists.at(unhandeledPt[0]) = dists.at(unhandeledPt[1]) + 1;
    handleQueue(ptQueue, basinLabels, gradients, dists, basinLabels(unhandeledPt[1]), gridSize, nDim);
  }

  // for( int i = 0; i < dists.size(); i++){
  //   if(i % gridSize == 0) Rcout << std::endl;
  //   Rcout << dists[i] << "";
  // }
  // Rcout << std::endl;
  return basinLabels;
}

