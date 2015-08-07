#ifndef __COUPLER_H__
#define __COUPLER_H__

#include "drift.h"
#include "omd/omdtool.h"
#include "omd/treader.h"

using namespace edrift;
using namespace omd;

// --------------- COUPLER CLASS ---------------

// Conductivity
class TeTa: public DriftCoupler {
  double k;
  double *te,*ta;

public:
  TeTa(double k0,double* tte,double* tta) {
    k=k0;
    te=tte;
    ta=tta;
  }

  void init(DriftSolver* ds) {
    DriftCoupler::init(ds);
    for(int ic=0;ic<solver->dim;ic++){
      if(solver->map[ic]) solver->con[ic]=k;
      else solver->con[ic]=0.0;
    }
  }

  void head() {
    double* con=solver->con;
    for(int i=0;i<solver->ncp;i++) {
      int idx=solver->clist[solver->me][i];
      if(!solver->map[idx]) continue; // skip unmapped cells
      con[idx]=ta[idx]>0?k*te[idx]/ta[idx]:0.0;
    }
  }

};

// Table reading, for example Heat-capacity
// how it works: for each cell, read table value at corresponding point stored in xptr.
// put the value into target_ptr. Example, if xptr is temperature, read a function of temperature
// of the table.
// use omd compact tablename syntax: table@filename:scale:constant

class TableFill: public DriftCoupler {
  double *a,*b;
  TableReader* table;

public:
  TableFill(string table_name,double* xptr, double* target_ptr){
    table=new TableReader(table_name);
    a=xptr; b=target_ptr;
  }

  ~TableFill(){delete table;}

  void head(){
    for(int i=0;i<solver->ncp;i++) {
      int idx=solver->clist[solver->me][i];
      if(!solver->map[idx]) continue; // skip unmapped cells
      b[idx]=table->read(a[idx]);
    }
  }

};

//* fill yptr as a linear function of xptr: yptr=factor*xptr
class LinearFill: public DriftCoupler {
  double *a,*b, fac;

public:
  LinearFill(double factor,double* xptr, double* target_ptr){
    fac=factor;
    a=xptr; b=target_ptr;
  }

  void head(){
    for(int i=0;i<solver->ncp;i++) {
      int idx=solver->clist[solver->me][i];
      if(!solver->map[idx]) continue; // skip unmapped cells
      b[idx]=fac*a[idx];
    }
  }

};

// electron phonon coupling

class TTMCoupling: public DriftCoupler {
  double G_const;
  double* Ta;
  double* Te;
  TableReader* G_table;

public:
  TTMCoupling(double coupling_factor, double* te, double* ta) {
    G_const=coupling_factor;
    Te=te;
    Ta=ta;
    G_table=NULL;
  }

  // coupling_table uses omd compact tablename syntax: table@filename:scale:constant
  TTMCoupling(string coupling_table, double* te, double* ta) {
    G_table=new TableReader(coupling_table);
    Ta=ta;
    Te=te;
  }

  ~TTMCoupling(){delete G_table;}

  void head(){
    for(int i=0;i<solver->ncp;i++) {
      int idx=solver->clist[solver->me][i];
      if(!solver->map[idx]) continue; // skip unmapped cells
      if(G_table)
        solver->Ext[idx]=G_table->read(Te[idx])*(Ta[idx]-Te[idx]);
      else
        solver->Ext[idx]=G_const*(Ta[idx]-Te[idx]);
    }
  }

};

// Propagating Lambert-Beer in +z direction
// TODO: another directions
// pulse function can be pulse_width of a normalized function from table.
// the source term is multiplied by the function.

class ProLB: public DriftCoupler {
  double A;
  double n0;
  double pw;
  double* dens;
  double* src;
  double del,lam;
  int nx,ny,nz;
  TableReader* ptab;

public:
  ProLB(double E_0, // fluence in energy/length^2
        double lambda,
        double equ_density,
        double pulse_width,
        double* density_map
        )
  {
    A=E_0*equ_density;
    lam=lambda;
    n0=equ_density;
    dens=density_map;
    pw=pulse_width;
    ptab=NULL;
  }

  ProLB(double E_0,
        double lambda,
        double equ_density,
        const char* pulse_table,  // in omd tablename format
        double* density_map
        )
  {
    // fluence A=E_0*n_0
    A=E_0*equ_density;
    lam=lambda;
    n0=equ_density;
    dens=density_map;
    ptab=new TableReader(pulse_table);
  }

  ~ProLB(){if(ptab) delete ptab;}

  void init(DriftSolver* ds) {
    DriftCoupler::init(ds);
    src=solver->GetSourceVector();
    nx=solver->nx;
    ny=solver->ny;
    nz=solver->nz;
  }

  void head(){ // for other direction only change this function
    double norm;
    double t=(double)solver->step * solver->timestep;

    if(ptab){ // table is only _time function_
      if(t>ptab->max_range()) return;
      norm=ptab->read(t);
    } else {
      if(t>pw) return;
      norm=(1./pw);
    }
    del=solver->GetSpacestep()/lam;

    for(int ix=0;ix<nx;ix++)
      for(int iy=0;iy<ny;iy++) {
        src[solver->index(ix,iy,0)]=A*norm;
        for(int iz=1;iz<nz;iz++) {
          int i0=solver->index(ix,iy,iz-1);
          int i1=solver->index(ix,iy,iz);
          double beta=dens[i1]/n0;
           src[i1]=src[i0]*exp(-beta*del);
        }
      }
  }

};

// Homogen source with table time-function

class HomogenSource: public DriftCoupler {
  double eps;
  double pw;
  double *src;
  TableReader *ptab;

public:

  // Strength is in energy density: Energy/volume
  HomogenSource(double strength, double pulse_width) {
    eps=strength; pw=pulse_width;ptab=NULL;
  }

  HomogenSource(double strength, string pulse_table) {
    eps=strength; ptab=new TableReader(pulse_table);
  }

  ~HomogenSource(){if(ptab) delete ptab;}

  void init(DriftSolver* ds) {
    DriftCoupler::init(ds);
    src=solver->GetSourceVector();
  }

  void head(){ // for other direction only change this function
    double eng;
    double t=(double)solver->step * solver->timestep;

    if(ptab){ // table is only _time function_
      if(t>ptab->max_range()) return;
      eng=eps*ptab->read(t);
    } else {
      if(t>pw) return;
      eng=(eps/pw);
    }

    for(int i=0;i<solver->ncp;i++) {
      int idx=solver->clist[solver->me][i];
      if(!solver->map[idx]) continue; // skip unmapped cells
      src[idx]=eng;
    }

  }
};

#endif
