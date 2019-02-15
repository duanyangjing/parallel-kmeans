#include <iostream>
#include <fstream>
#include <vector>
#include <pthread.h>
#include <math.h>


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
bool find_home_cluster(Point& p) {
  double mindist = p.dist(clusters[0]->centroid);
  Cluster* home = nullptr;
  for (Cluster* cluster : clusters) {
    double d = p.dist(cluster->centroid);
    if (d <= mindist) {
      mindist = d;
      home = cluster;
    }
  }
  home->members.push_back(&p);

  bool converge = p.i == home->i;
  p.i = home->i;

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
    *converge = *converge && find_home_cluster(*points[i]);
  }

  // free args before thread exits
  delete params;
  // equivalent to calling pthread_exit with return value
  pthread_exit((void *)converge);
}


// return if the assignment was converged (same as previous assignment)
bool assign() {
  int chunksize = points.size() / nthreads;
  for (int i = 0; i < nthreads; i++) {
    // dynamically allocate arguments, running thread should free it
    // before exit
    Args* args = new Args(i * chunksize, (i+1) * chunkSize);
    pthread_create(&threads[i], nullptr, assign_on_chunk, args);
  }

  bool converge = true;
  for (int i = 0; i < nthreads; i++) {
    void *ret;
    pthread_join(threads[i], &ret);
    converge = converge && (*(bool*)ret);
    // free return of a thread after accumulating its value
    delete ret;
  }

  return converge;
} 

struct Argss {
  int si, ei;
  double* sums;
  Argss(int si, int ei, double* sums):
    si(si), ei(ei), sums(sums) {}
};


void* update_centroid_on_chunk(void* argss) {
  Argss* params = (Argss*)argss;
  for (int i = params->si; i < params->ei && i < points.size(); i++) {
    Point p = points[i];
    int home = p.i;
    for (int d = 0; d < D; d++) {
      sums[home][d] += 
    }
  }
}

  
//
void update_centroid(int K) {
  // K points, each point has D dimension
  double* sums = new double[K];
  int chunksize = points.size() / nthreads;
  for (int i = 0; i < nthreads; i++) {
    Argss* argss = new Argss(i*chunksize, (i+1)*chunkSize, sums);
    pthread_create(&threads[i], nullptr, update_centroid_on_chunk, argss);
  }

  delete[] sums;
  delete[] counts;
}

void kmeans(int K) {
  for (int i = 0; i < K; i++) {
    Cluster* cluster = new Cluster(i, points[i]->coords);
    clusters.push_back(cluster);
  }
  bool converge = assign();
  int niters = 1;
  while(!converge && niters < 20) {
    
  }

}



int main(int argc, char** argv) {
  char* filename = argv[1];
  int K = atoi(argv[2]); // number of clusters
  nthreads = atoi(argv[3]);
  
  init(filename);
  kmeans(K);
  
  // for (Point* p : points) {
  //   for (auto x : p->coords) {
  //     std::cout << x << " ";
  //   }
  //   std::cout << "\n";
  // }
  
  return 0;
}
