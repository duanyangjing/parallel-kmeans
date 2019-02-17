#include "utils.hh"


Matrix::Matrix(int nr, int nc): nr(nr), nc(nc) {
  data = new float[nr * nc];
}
  
Matrix::~Matrix() {
  delete[] data;
}
  
float Matrix::get(int r, int c) {
  if (r * c >= nr * nc) throw "Matrix index out of bounds.";
  return data[r * nc + c];
}

float* Matrix::getRow(int r) {
  if (r < 0 || r >= nr) throw "Matrix index out of bounds.";
  return &data[r * nc];
}

void Matrix::set(int r, int c, float v) {
  if (r * c >= nr * nc) throw "Matrix index out of bounds.";
  data[r * nc + c] = v;
}

void Matrix::setRow(int r, float* p) {
  if (r < 0 || r >= nr) throw "Matrix index out of bounds.";
  for (int c = 0; c < nc; c++) {
    data[r * nc + c] = p[c];
  }
}
