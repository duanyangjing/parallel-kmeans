#include <pthread.h>
#include "utils.hh"

#define MAX_ITERATIONS 20

class ParallelKMeans;

struct Args {
  int si, ei; // si inclusive, ei exclusive
  ParallelKMeans* self;
  Args(int si, int ei, ParallelKMeans* self): si(si), ei(ei), self(self) {}
};


// struct Argss {
//   int si, ei;
//   Matrix* avgs;
//   ParallelKMeans* self;
//   Argss(int si, int ei, Matrix* avgs, ParallelKMeans* self):
//     si(si), ei(ei), avgs(avgs), self(self) {}
// };

class ParallelKMeans {
private:
  int N;  // number of points
  int D;  // dimensionality
  int K;  // number of clusters
  int nthreads;

  pthread_t* threads; // size nthreads array
  Matrix* points; // N * D matrix
  Matrix* clusters; // K * D matrix, centroids of clusters
  int* assignment; // size N array, map point p (row index of points) to cluster c (row index of clusters)
  int* population; // size K array, population of a cluster

  float dist(float* p1, float* p2);
  bool findHomeCluster(int pi);
  // This seems to be a workaround to call pthread on a class member function
  // more elegant solution should be using std::thread
  static void* assignChunkWrapper(void* args) {
    ParallelKMeans* obj = (ParallelKMeans*)(((Args*)args)->self);
    return obj->assignChunk(args);
  }
  void* assignChunk(void* args);
  bool assign();

  static void* updateCentroidsChunkWrapper(void* args) {
    ParallelKMeans* obj = (ParallelKMeans*)(((Args*)args)->self);
    return obj->updateCentroidsChunk(args);
  }
  void* updateCentroidsChunk(void* argss);
  void updateCentroids();
  

public:
  // stream for points
  ParallelKMeans(std::ifstream& in, int K, int nthreads);
  ~ParallelKMeans();
  void writeCentroids(std::ofstream& out);
  void writePointAssignment(std::ofstream& out);
  // run clustering
  void run();
};
