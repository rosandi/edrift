#include <omd/omdtool.h>
#include "value.h"

value::value(double* vptr, double v, const char* csym, int sta, int sto) {
  ptr=vptr;
	val=v;
	start=sta;
	stop=sto;
  lsym=strlen(csym);
  sprintf(sym,"%s",csym);	
}

void value::init(DriftSolver* ds) {
  DriftCoupler::init(ds);
}

void valuehead::head(){
  if(solver->step<start) return;
  if(solver->step>stop) return;
  
  for(int i=0;i<solver->dim;i++){
    char ch=solver->GetMapSymbol(i);
    int c=0;

    for(;c<lsym;c++) if(ch==sym[c]) break;
    if(c==lsym) continue;

    ptr[i]=val;

  }
}  
    //std::cout << "sym="<<ch<<" val=" << ptr[i] << "\n";

void valuetail::tail(){
  if(solver->step<start) return;
  if(solver->step>stop) return;
  
  for(int i=0;i<solver->dim;i++){
    char ch=solver->GetMapSymbol(i);
    int c=0;
    for(;c<lsym;c++) if(ch==sym[c]) break;
    if(c==lsym) continue;
    ptr[i]=val;
  }
}
