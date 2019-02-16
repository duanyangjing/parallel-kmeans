#include <iostream>
#include <fstream>
#include <vector>
#include <pthread.h>
#include <math.h>
#include "utils.h"
#include <assert.h>
#include <iomanip>


typedef std::vector<double> coords_t;

double coord_dist(const coords_t& p1, const coords_t& p2) {
  double sum = 0.0;
  for (int i = 0; i < p1.size(); i++) {
    sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);
  }
  
  return sqrt(sum);
}

struct Point {
  int i; // assignment to cluster id [0, K)
  coords_t coords;
  Point(coords_t coords): coords(coords) {i = -1;}
  double dist(coords_t coords) {
    return coord_dist(coords, this->coords);
  }
};

std::ostream& operator <<(std::ostream& os, const Point& p) {
  for (double coord : p.coords) {
    os << coord << " ";
  }
  return os;
}

typedef std::vector<Point*> points_t;

struct Cluster {
  int i; // cluster id [0, K) 
  coords_t centroid;
  points_t members;
  Cluster(int i, coords_t c): i(i), centroid(c) {}
};



int N, D; // N - number of points, D - dimensionality
int nthreads;
#define MAX_THREADS 16
#define MAX_ITERATIONS 20
pthread_t threads[MAX_THREADS];
std::vector<Cluster*> clusters;
points_t points;
pthread_mutex_t lock;


// read file and initialize parameters
void init(const char* file) {
  // initialize N, D, points
  std::ifstream in(file);
  in >> N >> D;
  
  for (int i = 0; i < N; i++) {
    coords_t* coords = new std::vector<double>();
    double x;
    for (int i = 0; i < D; i++) {
      in >> x;
      coords->push_back(x);
    }
    points.push_back(new Point(*coords));
  }
}


// find the closest cluster to given point p, add p to the members
// of the closest cluster. Concurrent read to centroid field, but
// write to members field, so no locks needed.
// return true if the point didn't change the cluster it belongs to
bool find_home_cluster(Point* p) {
  double mindist = p->dist(clusters[0]->centroid);
  Cluster* home = nullptr;
  for (Cluster* cluster : clusters) {
    double d = p->dist(cluster->centroid);
    if (d <= mindist) {
      mindist = d;
      home = cluster;
    }
  }
  home->members.push_back(p);

  bool converge = p->i == home->i;
  p->i = home->i;
  assert(p->i >= 0);

  pthread_mutex_lock(&lock);
  std::cout<<"Assigned point "<< *p << " to cluster " << p->i <<"\n";
  pthread_mutex_unlock(&lock);
  return converge;
}


struct Args {
  int si, ei; // si inclusive, ei exclusive
  Args(int si, int ei): si(si), ei(ei) {}
};


void* assign_on_chunk(void* args) {
  Args* params = (Args*)args;
  // dynamically allocate the return value, master thread frees it
  bool* converge = new bool;
  *converge = true;
  for (int i = params->si; i < params->ei && i < points.size(); i++) {
    pthread_mutex_lock(&lock);
    std::cout<<"Working on assigning point " << i <<"\n";
    pthread_mutex_unlock(&lock);
    // && is short circuit, put the function call at left side !!!
    *converge = find_home_cluster(points[i]) && *converge;
  }

  // free args before thread exits
  delete params;
  // equivalent to calling pthread_exit with return value
  pthread_exit((void *)converge);
}


// return if the assignment was converged (same as previous assignment)
bool assign() {
  // !! This must be ceiling, otherwise some points will be left out!!
  int chunksize = points.size() / nthreads + (points.size() % nthreads != 0);
  std::cout<<"Assign points to cluster... Divide work to chunksize "<<chunksize<<"\n";
  for (int i = 0; i < nthreads; i++) {
    int si = i * chunksize;
    int ei = (i+1) * chunksize;
    // dynamically allocate arguments, running thread should free it
    // before exit
    Args* args = new Args(si, ei);
    std::cout<<"generate work from "<< si << " to " << ei << "\n";
    pthread_create(&threads[i], nullptr, assign_on_chunk, args);
  }

  bool converge = true;
  for (int i = 0; i < nthreads; i++) {
    void *ret;
    pthread_join(threads[i], &ret);
    converge = converge && (*(bool*)ret);
    // free return of a thread after accumulating its value
    delete (bool*)ret;
  }

  std::cout<<"Assign done.\n";
  return converge;
} 

struct Argss {
  int si, ei;
  Matrix* avgs;
  Argss(int si, int ei, Matrix* avgs):
    si(si), ei(ei), avgs(avgs) {}
};


void* update_centroid_on_chunk(void* argss) {
  Argss* params = (Argss*)argss;
  for (int i = params->si; i < params->ei && i < points.size(); i++) {
    Point* p = points[i];
    assert(p->i >= 0);
    int home = p->i;
    std::cout<<"Updating centroid for cluster " << home << "\n";
    for (int d = 0; d < D; d++) {
      double oldavg = params->avgs->get((std::size_t)home, (std::size_t)d);
      params->avgs->set((std::size_t)home, (std::size_t)d, oldavg + p->coords[d] / clusters[home]->members.size()); 
    }
  }

  // return origin params, because matrix inside has been updated for return.
  pthread_exit((void*)params);
}

  
//
void update_centroid(int K) {
  std::cout<<"Update centroid...\n";
  // K clusters / points, each point has D dimension
  Matrix* avgs = new Matrix(K, D);
  for (int k = 0; k < K; k++) {
    for (int d = 0; d < D; d++) {
      avgs->set(k,d,0.0);
    }
  }
  int chunksize = points.size() / nthreads + (points.size() % nthreads != 0);
  for (int i = 0; i < nthreads; i++) {
    int si = i * chunksize;
    int ei = (i + 1) * chunksize;
    Argss* argss = new Argss(si, ei, avgs);
    std::cout<<"generate work from "<< si << " to " << ei << "\n";
    pthread_create(&threads[i], nullptr, update_centroid_on_chunk, argss);
  }

  // after join, matrix contains K new centroid
  for (int i = 0; i < nthreads; i++) {
    pthread_join(threads[i], nullptr);
  }

  // TODO: this can be parallized too
  for (int k = 0; k < K; k++) {
    for (int d = 0; d < D; d++) {
      (clusters[k]->centroid)[d] = avgs->get((std::size_t)k,(std::size_t)d);
    }
  }
  // TODO: possible mem leak, argss (si,ei) not freed, but very tiny
  delete avgs;
}

void kmeans(int K) {
  // initialize k clusters
  for (int i = 0; i < K; i++) {
    Cluster* cluster = new Cluster(i, points[i]->coords);
    clusters.push_back(cluster);
  }
  bool converge = assign();

  
  for (auto p : points) {
    assert(p->i >= 0);
  }
  
  
  int niters = 1;
  while(!converge && niters < 20) {
    update_centroid(K);
    niters++;
  }
}

void output(int K, const char* f1, const char* f2) {
  std::ofstream s1(f1);
  std::ofstream s2(f2);
  for (int i = 0; i < N; i++) {
    s1 << points[i]->i << "\n";
  }

  for (int k = 0; k < K; k++) {
    for (int d = 0; d < D; d++) {
      s2 << std::fixed << std::setprecision(4) << (clusters[k]->centroid)[d] << " ";
    }
    s2 << "\n";
  }
}


int main(int argc, char** argv) {
  char* filename = argv[1];
  int K = atoi(argv[2]); // number of clusters
  nthreads = atoi(argv[3]);
  
  init(filename);
  kmeans(K);

  output(K, "clusters.txt", "centroids.txt");
  
  // for (Point* p : points) {
  //   for (auto x : p->coords) {
  //     std::cout << x << " ";
  //   }
  //   std::cout << "\n";
  // }
  
  return 0;
}
