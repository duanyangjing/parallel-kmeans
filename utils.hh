#ifndef _UTIL_
#define _UTIL_

class Matrix {
private:
  float* data;
  int nr;
  int nc;

public:
  Matrix(int nr, int nc);
  
  ~Matrix();
    
  float get(int r, int c);

  float* getRow(int r);
 
  void set(int r, int c, float v);

  // copy given float array to row r in the matrix.
  void setRow(int r, float* p);
  
};


#endif
