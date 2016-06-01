#ifndef __VALUE_H__
#define __VALUE_H__

#include "drift.h"

using namespace edrift;

// set constant value 
class value: public DriftCoupler {
  friend double cell_value(double i, double j, double k);
  friend double src_value(double i, double j, double k);
  friend double ext_value(double i, double j, double k);

protected:
  char sym[64];
  int lsym;
  int start;
  int stop;
  double val;
  double *ptr;

public:
  value(double* vptr, double v, const char* csym, int sta, int sto);
  void init(DriftSolver* ds);
};

class valuehead: public value {
public:
  valuehead(double* vptr, double v,const char* csym, int sta, int sto):
  value(vptr,v,csym,sta,sto){}
  void head();
};

class valuetail: public value {
public:
  valuetail(double* vptr, double v,const char* csym, int sta, int sto):
  value(vptr,v,csym,sta,sto){}
  void tail();
};


#endif
