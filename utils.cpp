#include <cstddef>

class Matrix {
private:
  double* data;
  std::size_t nr;
  std::size_t nc;

public:
  Matrix(std::size_t nr, std::size_t nc): nr(nr), nc(nc) {
    data = new double[nr * nc];
  }
  
  ~Matrix() {
    delete[] data;
  }
  
  double get(std::size_t x, std::size_t y) {
    if (x * y >= nr * nc) throw "Matrix index out of bounds.";
    return data[x * nc + y];
  }
 
  void set(std::size_t x, std::size_t y, double v) {
    if (x * y >= nr * nc) throw "Matrix index out of bounds.";
    data[x * nc + y] = v;
  }
};
