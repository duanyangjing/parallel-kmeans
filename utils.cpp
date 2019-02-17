#include "utils.hh"


Matrix::Matrix(int nr, int nc): nr(nr), nc(nc) {
  data = new double[nr * nc];
}
  
Matrix::~Matrix() {
  delete[] data;
}
  
double Matrix::get(int r, int c) {
  if (r * c >= nr * nc) throw "Matrix index out of bounds.";
  return data[r * nc + c];
}

double* Matrix::getRow(int r) {
  if (r < 0 || r >= nr) throw "Matrix index out of bounds.";
  return &data[r * nc];
}

void Matrix::set(int r, int c, double v) {
  if (r * c >= nr * nc) throw "Matrix index out of bounds.";
  data[r * nc + c] = v;
}
