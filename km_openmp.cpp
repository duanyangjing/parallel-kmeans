#include <iostream>
#include <fstream>
#include <pthread.h>
#include <math.h>
#include <stdlib.h>
#include "utils.hh"
#include <assert.h>
#include <iomanip>
#include <omp.h>

#define MAX_ITERATIONS 20
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



//-----------------------------Global stuff-------------------------------------
int N;
int D;
int K;

int nthreads;
Matrix* points; // N * D matrix
Matrix* clusters; // K * D matrix, centroids of clusters  
  
int* assignment; // size N array, map point p (row index of points) to cluster c (row index of clusters)
int* population; // size K array, population of a cluster
omp_lock_t* locks; // size K array, protect each cluster when computing centroids. 

//-----------------------------Private methods----------------------------------
float dist(float* p1, float* p2) {
  float sum = 0.0;
  for (int i = 0; i < D; i++) {
    sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);
  }
  
  return sqrt(sum);
}

// return if the assignment was converged (same as previous assignment)
bool assign() {
  bool converge = true;
  #pragma omp parallel for reduction(&&:converge) num_threads(nthreads)  
  // find home cluster for every point, reduction on converge  
  for (int pi = 0; pi < N; pi++) {
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

    converge = converge && assignment[pi] == home;
    assignment[pi] = home;
    //assert(assignment[pi] >= 0 && home >= 0);
    #pragma omp atomic
    population[home]++;
  }
  
  return converge;
}

void updateCentroids() {
  // K clusters / centroids, each centroid has D dimension
  // zero fill clusters
  for (int k = 0; k < K; k++) {
    for (int d = 0; d < D; d++) {
      clusters->set(k, d, 0.0);
    }
  }

  Matrix* localClusters[nthreads];
  for (int i = 0; i < nthreads; i++) {
    localClusters[i] = new Matrix(K, D);
    for (int k = 0; k < K; k++) {
      for (int d = 0; d < D; d++) {
	localClusters[i]->set(k, d, 0.0);
      }
    }
  }
  
  #pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < N; i++) {
    for (int d = 0; d < D; d++) {
      float coord = points->get(i, d);
      int home = assignment[i];
      
      int tid = omp_get_thread_num();
      Matrix* localCluster = localClusters[tid];
      float old = localCluster->get(assignment[i], d);
      localCluster->set(assignment[i], d, old + coord);
    }
  }
  
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
void init(std::ifstream& in) {
  // initialize N, D, points, clusters, assignment
  in >> N >> D;
  points = new Matrix(N, D);
  clusters = new Matrix(K, D);
  assignment = new int[N];
  population = new int[K];
  locks = new omp_lock_t[K];

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
    omp_init_lock(&locks[i]);
  }
}

void cleanup() {
  delete points;
  delete clusters;  
  delete[] assignment;
  delete[] population;
  delete[] locks;
}

void writeCentroids(std::ofstream& out) {
  out << K << " " << D << "\n";
  for (int i = 0; i < K; i++) {
    for (int j = 0; j < D; j++) {
      out << std::fixed << std::setprecision(3) << clusters->get(i,j) << " ";
    }
    out << "\n";
  }
}

void writePointAssignment(std::ofstream& out) {
  for (int i = 0; i < N; i++) {
    out << assignment[i] << "\n";
  }
}


void run() {
  double ts = monotonic_seconds();

  assign();
  for (int i = 0; i < MAX_ITERATIONS; i++) {
    updateCentroids();
    if (assign()) break;
  }

  double te = monotonic_seconds();
  print_time(te - ts);
}


int main(int argc, char** argv) {
  char* f = argv[1];
  K = atoi(argv[2]); // number of clusters
  nthreads = atoi(argv[3]);

  std::ifstream in(f);
  init(in);

  run();
  
  std::ofstream out1("clusters.txt");
  std::ofstream out2("centroids.txt");
  writePointAssignment(out1);
  writeCentroids(out2);
  
  cleanup();
  return 0;
}
