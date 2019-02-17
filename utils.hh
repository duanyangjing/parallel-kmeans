#ifndef _UTIL_
#define _UTIL_

class Matrix {
private:
  double* data;
  int nr;
  int nc;

public:
  Matrix(int nr, int nc);
  
  ~Matrix();
    
  double get(int r, int c);

  double* getRow(int r);
 
  void set(int r, int c, double v);
  
};


#endif
