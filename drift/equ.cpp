#include <cstring>
#include <climits>
#include <omd/omdtool.h>
#include "equ.h"

using namespace omd;

DriftSolver *solverptr;
equ* equptr;

double cell_value(double i, double j, double k){
  int ii=(int)equptr->var_ix+(int)i;
  int jj=(int)equptr->var_iy+(int)j;
  int kk=(int)equptr->var_iz+(int)k;
  int idx=solverptr->index(ii,jj,kk);
  if(idx<0) {
    char err[512];
    sprintf(err,"out of range: (%d,%d,%d):%d",ii,jj,kk,idx);
    throw err;
  }
  return solverptr->A[idx];
}

double src_value(double i, double j, double k){
  int ii=(int)equptr->var_ix+(int)i;
  int jj=(int)equptr->var_iy+(int)j;
  int kk=(int)equptr->var_iz+(int)k;
  int idx=solverptr->index(ii,jj,kk);
  if(idx<0) {
    char err[512];
    sprintf(err,"out of range: (%d,%d,%d):%d",ii,jj,kk,idx);
    throw err;
  }
  return solverptr->S[idx];
}

double ext_value(double i, double j, double k){
  int ii=(int)equptr->var_ix+(int)i;
  int jj=(int)equptr->var_iy+(int)j;
  int kk=(int)equptr->var_iz+(int)k;
  int idx=solverptr->index(ii,jj,kk);
  if(idx<0) {
    char err[512];
    sprintf(err,"out of range: (%d,%d,%d):%d",ii,jj,kk,idx);
    throw err;
  }
  return solverptr->Ext[idx];
}

equ::equ(double* vptr, const char* sequ, const char* csym, int sta, int sto){
  ptr=vptr;

  istringstream sst(replace_char(remove_char(sequ,' '),";",' '));
  string ss;
  sst>>ss; // first token is the equation
  equpar.SetExpr(ss.c_str());

  vauxcnt=0;

  while(sst>>ss) {
    string vname;
    istringstream ssv(replace_char(ss,"=",' '));
    ssv>>vname>>vaux[vauxcnt++];
    equpar.DefineVar(vname.c_str(),&(vaux[vauxcnt++]));
    if(vauxcnt>=10) throw "ERROR@equ too many variables";
  }

  sprintf(sym,"%s",csym);
  lsym=strlen(csym);
  start=sta;
  stop=sto;
  equpar.DefineFun("RSOURCE", src_value); // read cell value at relative position
  equpar.DefineFun("RCELL", cell_value); // read cell value at relative position
  equpar.DefineFun("REXT", ext_value); // read cell value at relative position
  equpar.DefineConst("Pi",M_PI);
  equpar.DefineConst("kb",8.6173324E-5); // eV/K
  equpar.DefineVar("LX",&var_lx);
  equpar.DefineVar("LY",&var_ly);
  equpar.DefineVar("LZ",&var_lz);
  equpar.DefineVar("DT",&var_dt);
  equpar.DefineVar("DL",&var_dl);
  equpar.DefineVar("x",&var_x);
  equpar.DefineVar("y",&var_y);
  equpar.DefineVar("z",&var_z);
  equpar.DefineVar("t",&var_t);
  equpar.DefineVar("cell",&var_cell);
  equpar.DefineVar("src",&var_src);
  equpar.DefineVar("ext",&var_ext);
  equpar.DefineVar("ix",&var_ix);
  equpar.DefineVar("iy",&var_iy);
  equpar.DefineVar("iz",&var_iz);

}

void equ::init(DriftSolver* ds) {
  DriftCoupler::init(ds);
  update();
}

void equ::update() {
  var_lx=(double)solver->nx*solver->spacestep;
  var_ly=(double)solver->ny*solver->spacestep;
  var_lz=(double)solver->nz*solver->spacestep;
  var_dt=solver->timestep;
  var_dl=solver->spacestep;
}

void equhead::head(){
	std::cout<<"entry ("<<start<<")("<<stop<<")---";
  if(solver->step<start) return;
  if(solver->step>stop) return;
  std::cout<<"called\n";
  
  solverptr=solver;
  equptr=this;

  var_step=solver->step;
  var_t=var_step*var_dt;

  for(int i=0;i<solver->dim;i++){
    char ch=solver->GetMapSymbol(i);
    int c=0;

    for(;c<lsym;c++) if(ch==sym[c]) break;
    if(c==lsym) continue;

    try {
      var_cell=solver->A[i];
      var_src=solver->S[i];
      var_ext=solver->Ext[i];
      var_ix=double(i%solver->nx);
      var_iy=double(i/solver->nx);
      var_iz=double(i/(solver->nx*solver->ny));
      var_x=var_ix*var_dl;
      var_y=var_iy*var_dl;
      var_z=var_iz*var_dl;
      ptr[i]=equpar.Eval();
     //  std::cout << "/* "<<ptr[i] << std::endl;
    } catch (mu::Parser::exception_type &e) {
      std::cerr << e.GetMsg() << std::endl;
      throw "ERROR@equhead";
    }

  }
}

void equtail::tail(){
  if(solver->step<start) return;
  if(solver->step>stop) return;

  solverptr=solver;
  equptr=this;
  for(int i=0;i<solver->dim;i++){
    char ch=solver->GetMapSymbol(i);
    int c=0;

    for(;c<lsym;c++) if(ch==sym[c]) break;
    if(c==lsym) continue;

    try {
      var_cell=solver->A[i];
      var_src=solver->S[i];
      var_ext=solver->Ext[i];
      var_ix=double(i%solver->nx);
      var_iy=double(i/solver->nx);
      var_iz=double(i/(solver->nx*solver->ny));
      var_x=var_ix*var_dl;
      var_y=var_iy*var_dl;
      var_z=var_iz*var_dl;
      ptr[i]=equpar.Eval();
    } catch (mu::Parser::exception_type &e) {
      std::cerr << e.GetMsg() << std::endl;
      throw "ERROR@equhead";
    }

  }
}
