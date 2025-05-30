#include <iostream>
#include <fstream>
#include <pthread.h>
#include <math.h>
#include <stdlib.h>
#include "utils.hh"
#include <assert.h>
#include <iomanip>
#include "pkmeans.hh"


//--------------------------------------Timer-----------------------------------
/* Gives us high-resolution timers. */
#define _POSIX_C_SOURCE 199309L
#include <time.h>

/* OSX timer includes */
#ifdef __MACH__
  #include <mach/mach.h>
  #include <mach/mach_time.h>
#endif

/**
* @brief Return the number of seconds since an unspecified time (e.g., Unix
*        epoch). This is accomplished with a high-resolution monotonic timer,
*        suitable for performance timing.
*
* @return The number of seconds.
*/
static inline double monotonic_seconds()
{
#ifdef __MACH__
  /* OSX */
  static mach_timebase_info_data_t info;
  static double seconds_per_unit;
  if(seconds_per_unit == 0) {
    mach_timebase_info(&info);
    seconds_per_unit = (info.numer / info.denom) / 1e9;
  }
  return seconds_per_unit * mach_absolute_time();
#else
  /* Linux systems */
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
#endif
}

/**
* @brief Output the seconds elapsed while clustering.
*
* @param seconds Seconds spent on k-means clustering, excluding IO.
*/
static void print_time(double const seconds)
{
  printf("k-means clustering time: %0.04fs\n", seconds);
}


//-----------------------------Private methods----------------------------------
float ParallelKMeans::dist(float* p1, float* p2) {
  float sum = 0.0;
  for (int i = 0; i < D; i++) {
    sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);
  }
  
  return sqrt(sum);
}

// Find the closest cluster to given point index i. Concurrent read on clusters,
// but no write, Concurrent r/w on assignment, but is per point (thread).
// but concurrent incrementing population[home] which is not thread safe!!!
// return true if the point didn't change the cluster it belongs to
bool ParallelKMeans::findHomeCluster(int pi) {
  float* p = points->getRow(pi);
  float mindist = dist(p, clusters->getRow(0));
  int home = 0;
  for (int ci = 1; ci < K; ci++) {
    float* cluster = clusters->getRow(ci);
    float d = dist(p, cluster);
    if (d <= mindist) {
      mindist = d;
      home = ci;
    }
  }

  bool sameHome = assignment[pi] == home;
  assignment[pi] = home;
  //assert(assignment[pi] >= 0 && home >= 0);
  
  pthread_mutex_lock(&locks[home]);
  population[home]++;
  pthread_mutex_unlock(&locks[home]);

  return sameHome;
}


void* ParallelKMeans::assignChunk(void* args) {
  Args* params = (Args*)args;
  // dynamically allocate the return value, master thread frees it
  bool* converge = new bool;
  *converge = true;
  for (int pi = params->si; pi < params->ei && pi < N; pi++) {
    // && is short circuit, put the function call at left side !!!
    *converge = findHomeCluster(pi) && *converge;
  }

  // free args before thread exits
  delete params;
  // equivalent to calling pthread_exit with return value
  pthread_exit((void *)converge);
}


// return if the assignment was converged (same as previous assignment)
bool ParallelKMeans::assign() {
  // !! This must be ceiling, otherwise some points will be left out!!
  int chunksize = N / nthreads + (N % nthreads != 0);
  //std::cout<<"Assign points to cluster... Divide work to chunksize "<<chunksize<<"\n";
  for (int i = 0; i < nthreads; i++) {
    int si = i * chunksize;
    int ei = (i+1) * chunksize;
    // dynamically allocate arguments, which child thread should free before exit
    Args* args = new Args(si, ei, this);
    // std::cout<<"generate work from "<< si << " to " << ei << "\n";
    pthread_create(&threads[i], NULL, assignChunkWrapper, args);
  }

  bool converge = true;
  for (int i = 0; i < nthreads; i++) {
    void *ret;
    pthread_join(threads[i], &ret);
    converge = converge && (*(bool*)ret);
    // free return of a thread after accumulating its value
    delete (bool*)ret;
  }

  // std::cout<<"Assign done.\n";
  return converge;
}


// Concurrent r/w on clusters!!! not thread safe
void* ParallelKMeans::updateCentroidsChunk(void* args) {
  Args* params = (Args*)args;
  // local cluster of this chunk of work, to be returned
  Matrix* localCluster = new Matrix(K, D);
  for (int i = 0; i < K; i++) {
    for (int j = 0; j < D; j++) {
      localCluster->set(i, j, 0.0);
    }
  }
  
  for (int pi = params->si; pi < params->ei && pi < N; pi++) {
    float* p = points->getRow(pi);
    int home = assignment[pi];
    for (int d = 0; d < D; d++) {
      float old = localCluster->get(home, d);
      localCluster->set(home, d, old + p[d]);
    }
  }

  delete params;
  // return origin params, because matrix inside has been updated for return.
  // no need to explicitly return avgs, its still there
  pthread_exit((void*)localCluster);
}

void ParallelKMeans::updateCentroids() {
  // K clusters / points, each point has D dimension
  // zero fill clusters
  for (int k = 0; k < K; k++) {
    for (int d = 0; d < D; d++) {
      clusters->set(k, d, 0.0);
    }
  }
  
  int chunksize = N / nthreads + (N % nthreads != 0);
  for (int i = 0; i < nthreads; i++) {
    int si = i * chunksize;
    int ei = (i + 1) * chunksize;
    Args* args = new Args(si, ei, this);
    pthread_create(&threads[i], NULL, updateCentroidsChunkWrapper, args);
  }

  // after join, each thread return a local cluster, need to sum and compute avg
  Matrix* localClusters[nthreads];
  for (int i = 0; i < nthreads; i++) {
    void* localCluster;
    pthread_join(threads[i], &localCluster);
    localClusters[i] = (Matrix*)localCluster;
  }

  // TODO: this can be parallized too, might be able to just
  // use the original clusters, no need to copy around
  
  for (int k = 0; k < K; k++) {
    //assert(population[k] > 0);
    for (int d = 0; d < D; d++) {
      float sum = 0.0;
      for (int t = 0; t < nthreads; t++) {
	sum += localClusters[t]->get(k, d);
      }
      clusters->set(k, d, sum / population[k]);
    }
    // reset population
    population[k] = 0;
  }

  for (int i = 0; i < nthreads; i++) {
    delete localClusters[i];
  }
  
}


//-----------------------------Public methods-----------------------------------
ParallelKMeans::ParallelKMeans(std::ifstream& in, int K, int nthreads):
  K(K), nthreads(nthreads) {
  // initialize N, D, points, clusters, assignment
  in >> N >> D;
  threads = new pthread_t[nthreads];
  points = new Matrix(N, D);
  clusters = new Matrix(K, D);
  assignment = new int[N];
  population = new int[K];
  locks = new pthread_mutex_t[K];

  float x;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < D; j++) {
      in >> x;
      points->set(i, j, x);
    }
    assignment[i] = -1;
  }

  for (int i = 0; i < K; i++) {
    population[i] = 0;
    clusters->setRow(i, points->getRow(i));
    if (pthread_mutex_init(&locks[i], NULL) != 0) {
      printf("\n mutex init has failed\n");
    }
  }
}

ParallelKMeans::~ParallelKMeans() {
  delete points;
  delete clusters;
  delete[] threads;
  delete[] assignment;
  delete[] population;
  delete[] locks;
}

void ParallelKMeans::writeCentroids(std::ofstream& out) {
  out << K << " " << D << "\n";
  for (int i = 0; i < K; i++) {
    for (int j = 0; j < D; j++) {
      out << std::fixed << std::setprecision(3) << clusters->get(i,j) << " ";
    }
    out << "\n";
  }
}

void ParallelKMeans::writePointAssignment(std::ofstream& out) {
  for (int i = 0; i < N; i++) {
    out << assignment[i] << "\n";
  }
}


void ParallelKMeans::run() {
  double ts = monotonic_seconds();

  assign();
  for (int i = 0; i < MAX_ITERATIONS; i++) {
    updateCentroids();
    // stop if converges
    if (assign()) break;
  }

  double te = monotonic_seconds();
  print_time(te - ts);
}


int main(int argc, char** argv) {
  char* f = argv[1];
  int K = atoi(argv[2]); // number of clusters
  int nthreads = atoi(argv[3]);

  std::ifstream in(f);
  ParallelKMeans km(in, K, nthreads);

  km.run();
  
  std::ofstream out1("clusters.txt");
  std::ofstream out2("centroids.txt");
  km.writePointAssignment(out1);
  km.writeCentroids(out2);
  
  // for (Point* p : points) {
  //   for (auto x : p->coords) {
  //     std::cout << x << " ";
  //   }
  //   std::cout << "\n";
  // }
  
  return 0;
}
