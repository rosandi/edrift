#include <mpi.h>
#include <iostream>
#include "drift.h"
#include "coupler.h"
#include "equ.h"
#include "driftlib.h"

using namespace edrift;

void drift_create(int nx, int ny, int nz, char* sym, char* map, void** ptr) {
  int flag;
  MPI_Initialized(&flag);

  if (!flag) {
    int argc = 0;
    char **argv = NULL;
    MPI_Init(&argc,&argv);
  }
  int me=20, nproc=11;
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);

  DriftSolver* drf=new DriftSolver(nx,ny,nz,sym,map);
  *ptr=(void*) drf;
}

void drift_destroy(void* ptr) {
  DriftSolver *drf=(DriftSolver*)ptr;
  delete drf;
}

void drift_set_timestep(void* ptr,double timestep) {
  DriftSolver *drf=(DriftSolver*) ptr;
  drf->SetTimestep(timestep);
}

void drift_set_spacestep(void* ptr,double delta) {
  DriftSolver* drf=(DriftSolver*)ptr;
  drf->SetSpacestep(delta);
}

/* ------- ADD NEW ALGORITHMS HERE! ----------- */
void drift_set_algorithm(void* ptr,const char* algoname) {
  DriftSolver* drf=(DriftSolver*)ptr;
  if(strcmp(algoname,"basic")==0) {
    drf->SetAlgorithm(new DriftAlgorithm);
  } else {
    std::cerr<<"ERROR: Algorithm " << algoname << " not implemented"<<std::endl;
  }
}
/* --------------------------------------------*/

void drift_step(void* ptr,int nstep) {
  try {
    DriftSolver* drf=(DriftSolver*)ptr;
    drf->Step(nstep);
  } catch(const char* errmsg) {
    std::cerr<<"ERROR: "<<errmsg<<std::endl;
  }
}

void drift_dump(void* ptr,const char* fname) {
  DriftSolver* drf=(DriftSolver*)ptr;
  drf->Dump(fname);
}

void* drift_cell_vec(void* ptr) {
  DriftSolver* drf=(DriftSolver*)ptr;
  return (void*)drf->GetVector();
}

void* drift_source_vec(void* ptr) {
  DriftSolver* drf=(DriftSolver*)ptr;
  return (void*)drf->GetSourceVector();
}

void* drift_ext_vec(void* ptr) {
  DriftSolver* drf=(DriftSolver*)ptr;
  return (void*)drf->GetExtVector();
}

void* drift_con_vec(void* ptr) {
  DriftSolver* drf=(DriftSolver*)ptr;
  return (void*)drf->GetConductivity();
}

void* drift_cap_vec(void* ptr) {
  DriftSolver* drf=(DriftSolver*)ptr;
  return (void*)drf->GetCapacity();
}

void drift_setcell(void* ptr, double val, const char* st) {
  DriftSolver* drf=(DriftSolver*)ptr;
  drf->SetCell(st,val);
}

void drift_setsource(void* ptr, double val, const char* st) {
  DriftSolver* drf=(DriftSolver*)ptr;
  drf->SetSource(st,val);
}

void drift_setext(void* ptr, double val, const char* st) {
  DriftSolver* drf=(DriftSolver*)ptr;
  drf->SetExt(st,val);
}

void drift_setcon(void* ptr, double val, const char* st) {
  DriftSolver* drf=(DriftSolver*)ptr;
  drf->SetConductivity(st,val);
}

void drift_setcap(void* ptr, double val, const char* st) {
  DriftSolver* drf=(DriftSolver*)ptr;
  drf->SetCapacity(st,val);
}

void drift_collect(void* ptr, double* arr) {
  DriftSolver* drf=(DriftSolver*)ptr;
  drf->collect(arr);
}

double drift_getdim(void* ptr,int axis) {
  DriftSolver* drf=(DriftSolver*)ptr;
  int retv=0;
  if(axis==0) retv=drf->nx;
  if(axis==1) retv=drf->ny;
  if(axis==2) retv=drf->nz;
  return retv;
}

/* -------- COUPLER LIST ---------- */
void drift_equ_head(void* ptr, double* data, const char* equ, const char* cells, int sta, int sto) {
  DriftSolver* drf=(DriftSolver*)ptr;
  drf->LinkCoupler(new equhead(data,equ,cells,sta,sto));
}

void drift_equ_tail(void* ptr, double* data, const char* equ, const char* cells, int sta, int sto) {
  DriftSolver* drf=(DriftSolver*)ptr;
  drf->LinkCoupler(new equtail(data,equ,cells,sta,sto));
}

void drift_value_head(void* ptr, double* data, double val, const char* cells, int sta, int sto) {
  DriftSolver* drf=(DriftSolver*)ptr;
  drf->LinkCoupler(new valuehead(data,val,cells,sta,sto));	
}

void drift_value_tail(void* ptr, double* data, double val, const char* cells, int sta, int sto) {
  DriftSolver* drf=(DriftSolver*)ptr;
  drf->LinkCoupler(new valuetail(data,val,cells,sta,sto));	
}

void drift_teta(void* ptr, double k, double* te, double* ta) {
  DriftSolver* drf=(DriftSolver*)ptr;
  drf->LinkCoupler(new TeTa(k,te,ta));
}

void drift_table(void* ptr, const char* tname, double* x, double* y) {
  DriftSolver* drf=(DriftSolver*)ptr;
  drf->LinkCoupler(new TableFill(tname,x,y));
}

void drift_ttm_coupling_table(void* ptr, const char* tname, double* te, double* ta) {
  DriftSolver* drf=(DriftSolver*)ptr;
  drf->LinkCoupler(new TTMCoupling(tname,te,ta));
}

void drift_ttm_coupling(void* ptr, double cval, double* te, double* ta) {
  DriftSolver* drf=(DriftSolver*)ptr;
  drf->LinkCoupler(new TTMCoupling(cval,te,ta));
}

void drift_prolb(void* ptr, double e0, double lamb, double n0,  double* nmap, double tau) {
  DriftSolver* drf=(DriftSolver*)ptr;
  drf->LinkCoupler(new ProLB(e0,lamb,n0,tau,nmap));
}

void drift_prolb_table(void* ptr, double e0, double lamb, double n0, double* nmap, const char* tname) {
  DriftSolver* drf=(DriftSolver*)ptr;
  drf->LinkCoupler(new ProLB(e0,lamb,n0,tname,nmap));
}

void drift_source_pulse(void* ptr, double e0,double pwidth) {
  DriftSolver* drf=(DriftSolver*)ptr;
  drf->LinkCoupler(new HomogenSource(e0,pwidth));
}

void drift_source_table(void* ptr, double e0, const char* tname) {
  DriftSolver* drf=(DriftSolver*)ptr;
  drf->LinkCoupler(new HomogenSource(e0,tname));
}
