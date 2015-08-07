#ifndef __EQU_H__
#define __EQU_H__

#include <muParser/muParser.h>
#include "drift.h"

using namespace edrift;

class equ: public DriftCoupler {
  friend double cell_value(double i, double j, double k);
  friend double src_value(double i, double j, double k);
  friend double ext_value(double i, double j, double k);

protected:
  char sym[64];
  int lsym;
  double* ptr;
  int start;
  int stop;

  mu::Parser equpar;

  // equation variables
  double var_x;
  double var_y;
  double var_z;
  double var_t;
  double var_lx;
  double var_ly;
  double var_lz;
  double var_dt;
  double var_dl;
  double var_cell;
  double var_src;
  double var_ext;
  double var_step;
  double var_ix;
  double var_iy;
  double var_iz;
  double vaux[10];
  int vauxcnt;

public:
  equ(double* vptr, const char* sequ, const char* csym, int sta, int sto);
  void init(DriftSolver* ds);
  void update();
};

class equhead: public equ {
public:
  equhead(double* vptr, const char* sequ,const char* csym, int sta, int sto):
  equ(vptr,sequ,csym,sta,sto){}
  void head();
};

class equtail: public equ {
public:
  equtail(double* vptr, const char* sequ,const char* csym, int sta, int sto):
  equ(vptr,sequ,csym,sta,sto){}
  void tail();
};

#endif
