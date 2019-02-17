#include <iostream>
#include <fstream>
#include <pthread.h>
#include <math.h>
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
double ParallelKMeans::dist(double* p1, double* p2) {
  double sum = 0.0;
  for (int i = 0; i < D; i++) {
    sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);
  }
  
  return sqrt(sum);
}

// Find the closest cluster to given point index i. Concurrent read on clusters,
// but no write, Concurrent r/w on assignment, but is per point (thread).
// So no lock needed.
// return true if the point didn't change the cluster it belongs to
bool ParallelKMeans::findHomeCluster(int pi) {
  double* p = points->getRow(pi);
  double mindist = dist(p, clusters->getRow(0));
  int home = -1;
  for (int ci = 0; ci < K; ci++) {
    double* cluster = clusters->getRow(ci);
    double d = dist(p, cluster);
    if (d <= mindist) {
      mindist = d;
      home = ci;
    }
  }

  bool converge = assignment[pi] == home;
  assignment[pi] = home;
  assert(assignment[pi] >= 0 && home >= 0);
  population[home]++;

  return converge;
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
    //std::cout<<"generate work from "<< si << " to " << ei << "\n";
    pthread_create(&threads[i], nullptr, assignChunkWrapper, args);
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



void* ParallelKMeans::updateCentroidsChunk(void* argss) {
  Argss* params = (Argss*)argss;
  for (int pi = params->si; pi < params->ei && pi < N; pi++) {
    double* p = points->getRow(pi);
    int home = assignment[pi];
    for (int d = 0; d < D; d++) {
      double oldavg = params->avgs->get(home, d);
      params->avgs->set(home, d, oldavg + p[d] / population[home]); 
    }
  }

  // return origin params, because matrix inside has been updated for return.
  // no need to explicitly return avgs, its still there
  pthread_exit(nullptr);
}

void ParallelKMeans::updateCentroids() {
  // K clusters / points, each point has D dimension
  Matrix* avgs = new Matrix(K, D);
  for (int k = 0; k < K; k++) {
    for (int d = 0; d < D; d++) {
      avgs->set(k,d,0.0);
    }
  }
  
  int chunksize = N / nthreads + (N % nthreads != 0);
  for (int i = 0; i < nthreads; i++) {
    int si = i * chunksize;
    int ei = (i + 1) * chunksize;
    Argss* argss = new Argss(si, ei, avgs, this);
    pthread_create(&threads[i], nullptr, updateCentroidsChunkWrapper, argss);
  }

  // after join, matrix contains K new centroid
  for (int i = 0; i < nthreads; i++) {
    pthread_join(threads[i], nullptr);
  }

  // TODO: this can be parallized too, might be able to just
  // use the original clusters, no need to copy around
  for (int k = 0; k < K; k++) {
    for (int d = 0; d < D; d++) {
      clusters->set(k, d, avgs->get(k, d));
    }
  }
  // TODO: possible mem leak, argss (si,ei) not freed, but very tiny
  delete avgs;
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

  for (int i = 0; i < N; i++) {
    assignment[i] = -1;
  }

  double x;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < D; j++) {
      in >> x;
      points->set(i, j, x);
    }
  }
}

ParallelKMeans::~ParallelKMeans() {
  delete points;
  delete clusters;
  delete[] assignment;
  delete[] population;
}

void ParallelKMeans::writeCentroids(std::ofstream& out) {
  for (int i = 0; i < K; i++) {
    for (int j = 0; j < D; j++) {
      out << clusters->get(i,j) << " ";
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
  bool converge = assign();
  int niters = 1;

  double ts = monotonic_seconds();

  while (!converge && niters < 20) {
    updateCentroids();
    niters++;
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
