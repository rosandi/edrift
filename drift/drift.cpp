#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <climits>
#include <cfloat>
#include <mpi.h>
#include "omd/omdtool.h"
#include "drift.h"
// For 2D problems...

using std::vector;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::cout;
using std::endl;

using namespace omd;
using namespace edrift;
// equation parser
DriftSolver* myptr;
static char err_msg[1024];

//----------------------------------------------------------------
// on return str will be cleared if out of time range or no range

DriftSolver::DriftSolver(int _nx, int _ny, int _nz, string _symbol, string _map) {
  nx=_nx;
  ny=_ny;
  nz=_nz; // if nz<=2 then 2D simulation
  if(nz<=2) nz=1;

  dim=nx*ny*nz;
  strict_courant=true;
  memset(err_msg,0,1024);
  algo=NULL;
  map=NULL;
  bound=NULL;
  T=A=S=Ext=NULL;
  dE=dS=dExt=sumE=sumS=NULL;
  cap=con=NULL;
  ncoupler=0;
  step=0;
  timestep=1.0;
  spacestep=1.0;
  name="drift-solver";
  link_updated=true;
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  initiate_mem();
  Remap(_symbol.c_str(),_map.c_str());
}

void DriftSolver::initiate_mem() {

  // WARNING! dE is initiated on demand
  T=new double[dim];
  A=new double[dim];
  S=new double[dim];
  Ext=new double[dim];
  dS=new double[dim];
  dExt=new double[dim];
  con=new double[dim];
  cap=new double[dim];
  map = new int[dim];

  for(int i=0;i<dim;i++) {
    T[i]=A[i]=S[i]=Ext[i]=0.0;
    con[i]=cap[i]=1.0;
    map[i]=0; // all empty@init
  }

  counts=new int[nproc];
  displ=new int[nproc];
  recvnum=new int[nproc];
  sendnum=new int[nproc];
  lim=new int[2*nproc]; // lo&hi limits
  clist=new int*[nproc];
  nlist=new int*[nproc];
  glist=new int*[nproc];
  ndata=new double*[nproc];
  gdata=new double*[nproc];

  for(int p=0;p<nproc;p++) {
    clist[p]=NULL;
    nlist[p]=NULL;
    glist[p]=NULL;
    ndata[p]=NULL;
    gdata[p]=NULL;
  }

}

DriftSolver::~DriftSolver() {
  if(map) delete[] map;
  if(T) delete[] T;
  if(A) delete[] A;
  if(S) delete[] S;
  if(Ext) delete[] Ext;
  if(dE) delete[] dE;
  if(dS) delete[] dS;
  if(dExt) delete[] dExt;
  if(sumE) delete[] sumE;
  if(sumS) delete[] sumS;
  if(cap) delete[] cap;
  if(con) delete[] con;
  if(bound) delete[] bound;

  delete[] counts;
  delete[] displ;
  delete[] sendnum;
  delete[] recvnum;
  delete[] lim;

  if(nproc>1) {
    for (int p=0;p<nproc;p++) {
      if(ndata[p]){delete[] ndata[p];ndata[p]=NULL;}
      if(gdata[p]){delete[] gdata[p];gdata[p]=NULL;}
      if(nlist[p]){delete[] nlist[p];nlist[p]=NULL;}
      if(glist[p]){delete[] glist[p];glist[p]=NULL;}
      if(clist[p]){delete[] clist[p];clist[p]=NULL;}
    }

    delete[] ndata;
    delete[] gdata;
    delete[] nlist;
    delete[] glist;
    delete[] clist;
  }
}

bool DriftSolver::mycell(int ix, int iy, int iz) {
  int idx=index(ix,iy,iz);
  if(idx<0) return false;
  return (idx>=lim[me*2] && idx<=lim[me*2+1]);
}

void DriftSolver::select_ncell(int cond, int idx, char** nmap) {
  if(cond) {
    for(int p=0;p<nproc;p++) {
      if(p==me) continue;
      if(idx>=lim[p*2] && idx<=lim[p*2+1]) {
        nmap[p][idx]=1;
      }
    }
  }
}

void DriftSolver::send_receive_list() {
  MPI_Request req[2*nproc];
  MPI_Status  stat[2*nproc];
  int rq=0;

  for(int p=0;p<nproc;p++) { // send to p
    if(p==me) continue;
    if(recvnum[p]) {
      MPI_Send_init(glist[p], recvnum[p], MPI_INT, p, 1, MPI_COMM_WORLD, &req[rq++]);
    }
  }

  // fetch index list to send
  for(int p=0;p<nproc;p++) { // receive from p
    if(p==me) continue;
    if(sendnum[p]) {
      MPI_Recv_init(nlist[p], sendnum[p], MPI_INT, p, 1, MPI_COMM_WORLD, &req[rq++]);
    }
  }

  MPI_Startall(rq,req);
  MPI_Waitall(rq,req,stat);
  for(int i=0;i<rq;i++) MPI_Request_free(&req[i]);
}

// assign -1 on neigh empty
// empty neighbour subject to von neumann condition
// in this case the temperature is isolated in the target

void DriftSolver::get_neigcell_index(int idx, int& l, int& r, int& t, int& b, int& n, int& f) {
  int i,j,k;

  get_cellcoord(idx,i,j,k);
  CellBoundary *cc=&(bound[map[idx]]);
  l=i-1; r=i+1; t=j-1; b=j+1; n=k-1; f=k+1;

  if(l<0) {
    if(forced_periodic&1) l=nx-1;
    else if(cc->l==FDCELL_PERIODIC) l=nx-1;
  }

  if(r>=nx) {
    if(forced_periodic&1) r=0;
    else if(cc->r==FDCELL_PERIODIC) r=0;
  }
  if(t<0) {
    if(forced_periodic&2) t=ny-1;
    else if(cc->t==FDCELL_PERIODIC) t=ny-1;
  }
  if(b>=ny) {
    if(forced_periodic&2) b=0;
    else if(cc->b==FDCELL_PERIODIC) b=0;
  }

  if(nz!=1) {
    if(n<0) {
      if(forced_periodic&4) n=nz-1;
      else if(cc->n==FDCELL_PERIODIC) n=nz-1;
    }

    if(f>=nz){
      if(forced_periodic&4) f=0;
      else if(cc->f==FDCELL_PERIODIC) f=0;
    }
  } else n=f=0;

  l=cc->l?index(l,j,k):-1;
  r=cc->r?index(r,j,k):-1;
  t=cc->t?index(i,t,k):-1;
  b=cc->b?index(i,b,k):-1;

  if(nz!=1) {
    n=cc->n?index(i,j,n):-1;
    f=cc->f?index(i,j,f):-1;
  }

}

void DriftSolver::initiate_comm() {

  displ[0]=0;
  occup=0; // # cells pre proc
  for(int i=0;i<dim;i++) if(map[i])occup++;
  ncp=occup/nproc;

  if(me==nproc-1) { // last proc take responsible for tail
    ncp=occup-ncp*(nproc-1);
  }

  MPI_Allgather(&ncp,1,MPI_INT,counts,1,MPI_INT,MPI_COMM_WORLD);

  for(int i=1;i<nproc;i++) displ[i]=counts[i-1]+displ[i-1];
  for(int i=0;i<nproc;i++) clist[i]=new int[counts[i]];

  int hilo[2];

  {
    int nc=0;
    int cc=0;
    int i;
    int start=0;
    for(i=0;i<me;i++) start+=counts[i];

    for(i=0;i<dim;i++) {
      if(map[i]) {
        if((cc++)<start) continue;
        if(!nc) hilo[0]=i;
        clist[me][nc++]=i;
      }
      if(nc==ncp) break;
    }
    hilo[1]=i;
  }

  if(nproc==1) {
    lim[0]=hilo[0];
    lim[1]=hilo[1];
  }

  if(nproc>1) {

    {
      int tmp[occup];
      MPI_Allgatherv(clist[me],ncp,MPI_INT,tmp,counts,displ,MPI_INT,MPI_COMM_WORLD);

      int ii=0;
      for(int p=0;p<nproc;p++) {
        if(p==me) {ii+=ncp;continue;}
        for(int i=0;i<counts[p];i++) clist[p][i]=tmp[ii++];
      }
    }

    MPI_Allgather(hilo,2,MPI_INT,lim,2,MPI_INT,MPI_COMM_WORLD);

    char** nmap;
    nmap=new char*[nproc];
    for(int p=0;p<nproc;p++) {
      nmap[p]=new char[dim];
      for(int i=0;i<dim;i++) nmap[p][i]=0;
    }

    // check neighbor owner
    // avoid multiple-include
    for(int i=0;i<ncp;i++) {
      int idx=clist[me][i];
      int l,r,t,b,n,f;
//      CellBoundary *cc=&(bound[map[idx]]);
      get_neigcell_index(idx,l,r,t,b,n,f);

      select_ncell(l>=0,l,nmap);
      select_ncell(r>=0,r,nmap);
      select_ncell(t>=0,t,nmap);
      select_ncell(b>=0,b,nmap);

      if(nz!=1) {
        select_ncell(n>=0,n,nmap);
        select_ncell(f>=0,f,nmap);
      }

    }

    for(int p=0;p<nproc;p++) {
      recvnum[p]=0;
      sendnum[p]=0;
      for(int i=0;i<dim;i++)recvnum[p]+=(int)nmap[p][i];
    }

    for(int p=0;p<nproc;p++) {
      if(recvnum[p]==0||p==me) continue;
      glist[p]=new int[recvnum[p]];
      gdata[p]=new double[recvnum[p]];
      int n=0;
      for(int i=0;i<dim;i++)
        if(nmap[p][i]) glist[p][n++]=i;
    }

    // free...
    for(int p=0;p<nproc;p++) delete nmap[p];
    delete[] nmap;

    // ready and sync

    int allnum[nproc*nproc];
    MPI_Allgather(recvnum,nproc, MPI_INT, allnum, nproc, MPI_INT, MPI_COMM_WORLD);

    for(int p=0;p<nproc;p++) {
      if(p==me) continue;
      int* ptr=allnum+(p*nproc);
      sendnum[p]=ptr[me];
    }

    // allocate send and index buffer
    for(int p=0;p<nproc;p++) {
      if(!sendnum[p]) continue;
      nlist[p]=new int[sendnum[p]];
      ndata[p]=new double[sendnum[p]];
    }

    send_receive_list();
  }
}

// update couplers
void DriftSolver::update() {
  for(int nc=0;nc<ncoupler;nc++) coupler[nc]->update();
}

// remaping must be accompanied by setting/rereading cap&con
// if needed also src
// if the non-time-ranged cap/con expression, this sould be
// done automatic (sync_data).

void DriftSolver::Remap(const char* newsym, const char* newmap) {

  if(nproc>1) {
    for (int p=0;p<nproc;p++) {
      if(ndata[p]){delete[] ndata[p];ndata[p]=NULL;}
      if(gdata[p]){delete[] gdata[p];gdata[p]=NULL;}
      if(nlist[p]){delete[] nlist[p];nlist[p]=NULL;}
      if(glist[p]){delete[] glist[p];glist[p]=NULL;}
      if(clist[p]){delete[] clist[p];clist[p]=NULL;}
    }
  }

  // empty newsym or newmap -> leave untouched...
  if(newsym) {
    if(bound) delete[] bound;
    SetSymbolMap(newsym);
  }
  // FIXME!
  // WARNING! initiate only if map is changed...
  if(newmap) {
    SetMap(newmap);
    initiate_comm();
  }

}

void DriftSolver::SetTimestep(double dt){
  timestep=dt;
}

void DriftSolver::SetSpacestep(double h){
  spacestep=h;
  h2inv=1.0/(h*h);
}

int DriftSolver::search_symbol(char c) {
  for(int i=0;i<nbound;i++) {
    if(bound[i].sym==c) return i;
  }
  return 0; // @head -> empty cell '.'
}

char DriftSolver::map_symbol(int num) {
  if(num<0) return '.';
  if(num>=nbound) return '.';
  return bound[num].sym;
}

// example:
//  a ++++++

void DriftSolver::SetSymbolMap(string str) {
  istringstream ss(str.c_str());
  vector<CellBoundary> vb;
  CellBoundary cc;

  cc.sym='.'; // dot is used for empty cell: head of symbol list
  cc.l=cc.r=cc.t=cc.b=cc.n=cc.f=FDCELL_EMPTY;
  vb.push_back(cc);

  while(ss.good()) {
    string tok, sb;

    if(ss>>tok>>sb) {
      if(sb.size()<4)
        throw "invalid map symbol entry";

      for(int i=0;i<(int)tok.size();i++) {
        cc.sym=tok[i];
        cc.l=cc.r=cc.t=cc.b=cc.n=cc.f=FDCELL_EMPTY;

        if(sb[0]=='+') cc.l=FDCELL_OCCUPIED;
        if(sb[0]=='p') cc.l=FDCELL_PERIODIC;
        if(sb[1]=='+') cc.r=FDCELL_OCCUPIED;
        if(sb[1]=='p') cc.r=FDCELL_PERIODIC;
        if(sb[2]=='+') cc.t=FDCELL_OCCUPIED;
        if(sb[2]=='p') cc.t=FDCELL_PERIODIC;
        if(sb[3]=='+') cc.b=FDCELL_OCCUPIED;
        if(sb[3]=='p') cc.b=FDCELL_PERIODIC;
        if(nz>1 && sb.size()==6) { // 3D
          if(sb[4]=='+') cc.n=FDCELL_OCCUPIED;
          if(sb[4]=='p') cc.n=FDCELL_PERIODIC;
          if(sb[5]=='+') cc.f=FDCELL_OCCUPIED;
          if(sb[5]=='p') cc.f=FDCELL_PERIODIC;
        }

        vb.push_back(cc);
      }
    }
  }

  nbound=(int)vb.size();
  bound=new CellBoundary[nbound];
  for(int i=0;i<nbound;i++) {
    bound[i]=vb[i];
  }
/* DEBUG  
  for(int i=0;i<nbound;i++) {
	  std::cout<<bound[i].sym<<bound[i].l<<bound[i].r<<bound[i].t<<bound[i].b<<bound[i].n<<bound[i].f<<std::endl;
  }
*/

}

void DriftSolver::SetMap(const char* str) {
  for(int i=0;i<dim;i++) {
    if(str[i]==0x0) throw "insufficient number of character in map string";
    map[i]=search_symbol(str[i]);
  }
}

void DriftSolver::EnableDelta() {
  if(dE) delete[] dE;
  dE=new double[dim];
}

void DriftSolver::EnableTally() {
  if(sumE) delete[] sumE;
  if(sumS) delete[] sumS;
  sumE=new double[dim];
  sumS=new double[dim];
  memset(sumE,0,dim*sizeof(double));
  memset(sumS,0,dim*sizeof(double));
}

void DriftSolver::KeepSource() {
  if(dS) delete[] dS;
  if(dExt) delete[] dExt;
  dS=new double[dim];
  dExt=new double[dim];
  memset(dS,0,dim*sizeof(double));
  memset(dExt,0,dim*sizeof(double));
}

/**
 communicates data at neighboring borders
 */

void DriftSolver::communicate(double* ptr) {
  if(nproc>1) {

    MPI_Request req[2*nproc];
    MPI_Status  stat[2*nproc];
    int rq=0;

    for(int p=0;p<nproc;p++) {
      if(!sendnum[p]) continue;
      for(int i=0;i<sendnum[p];i++) ndata[p][i]=ptr[nlist[p][i]];
      MPI_Send_init(ndata[p],sendnum[p],MPI_DOUBLE,p,1,MPI_COMM_WORLD,&req[rq++]);
    }

    for(int p=0;p<nproc;p++) {
      if(recvnum[p])
        MPI_Recv_init(gdata[p],recvnum[p],MPI_DOUBLE,p,1,MPI_COMM_WORLD,&req[rq++]);
    }

    MPI_Startall(rq,req);
    MPI_Waitall(rq,req,stat);

    for(int p=0;p<nproc;p++) {
      if(recvnum[p]) for(int i=0;i<recvnum[p];i++) ptr[glist[p][i]]=gdata[p][i];
    }

    for(int i=0;i<rq;i++) MPI_Request_free(&req[i]);

  }
}

/**
 collects data and stor in every processor's array
 */

void DriftSolver::collect(double* ptr) {
  if(nproc>1) {
    double *tmp=new double[occup];
    double *dat=new double[ncp];

    for(int i=0;i<ncp;i++) dat[i]=ptr[clist[me][i]];

    MPI_Allgatherv(dat,ncp,MPI_DOUBLE,tmp,counts,displ,MPI_DOUBLE,MPI_COMM_WORLD);

    int ii=0;
    for(int p=0;p<nproc;p++) {
      if(p==me){ii+=ncp;continue;}
      for(int i=0;i<counts[p];i++) ptr[clist[p][i]]=tmp[ii++];
    }

    delete[] tmp;
    delete[] dat;
  }
}

void DriftSolver::collect() {
  collect(A);
  collect(dS);
  collect(dExt);
  if(dE) collect(dE);
  if(sumE) {
    collect(sumE);
    collect(sumS);
  }
}

void DriftSolver::GetIndexList(const char* sym, vector<int>& vind) {
  vind.clear();
  int ln=(int)strlen(sym);
  for(int i=0;i<dim;i++) {
    char ch=GetMapSymbol(i);
    for(int c=0;c<ln;c++) if(ch==sym[c]) vind.push_back(i);
  }
}

void DriftSolver::AssignArray(double* ptr, int ix, int iy, int iz, double val) {
  int idx=index(ix,iy,iz);
  if(idx<0) return;
  ptr[idx]=val;
}

void DriftSolver::AssignArray(double* ptr, const char* str, double val) {
  int ln=strlen(str);
  for(int i=0;i<dim;i++){
    char ch=GetMapSymbol(i);
    for(int c=0;c<ln;c++) if(ch==str[c]) ptr[i]=val;
  }
}

void DriftSolver::sync_data() {
  communicate(con);
  communicate(cap);
  communicate(A);
}

void DriftSolver::Step(int n) {
  // these are accumulated per call values
  memset(dS,0,dim*sizeof(double));
  memset(dExt,0,dim*sizeof(double));

  for(int loop=0;loop<n;loop++) {
    sync_data();
    for(int nc=0;nc<ncoupler;nc++) coupler[nc]->head();

    if(algo) {
      memcpy(T,A,sizeof(double)*dim);
      algo->exec(); // NULL means no diffusion
    }

    // factor 2 is from averaging of conductivity
    //
    // Ext array is used for coupling purpose.

    for(int i=0;i<ncp;i++) {

      int idx=clist[me][i];
      if(!map[idx]) continue; // skip unmapped cells

      if(cap[idx]<=0.0) {
        string msg("Heat capacity may not be zero at this point (drift.cpp@1228)!! ");
        msg.append(as_string(idx)+" "+as_string(map[idx])+" "+as_string(cap[idx]));
        throw msg.c_str();
      }

      if(algo) {
        if(dE) {
          dE[idx]=(2.0*A[idx]*h2inv+S[idx]+Ext[idx])*timestep; // gives dE/vol
          A[idx]=T[idx]+dE[idx]/cap[idx];
        } else {
          A[idx]=T[idx]+(2.0*A[idx]*h2inv+S[idx]+Ext[idx])*timestep/cap[idx];
        }
      } else {
        if(dE) {
          dE[idx]=(S[idx]+Ext[idx])*timestep; // gives dE/vol
          A[idx]+=dE[idx]/cap[idx];
        } else {
          A[idx]+=(S[idx]+Ext[idx])*timestep/cap[idx];
        }
      }

      // in eV/A^3
      dS[idx]+=S[idx]*timestep;
      dExt[idx]+=Ext[idx]*timestep;

    }

    // reserved capability to call a coupler class


    step++;
  } /* for loop */

  if(sumE) { // tally: if enabled
    for(int i=0;i<ncp;i++) {
      int idx=clist[me][i];
      sumS[idx]+=dS[idx];
      sumE[idx]+=dExt[idx];
    }
  }
}

double DriftSolver::GetSum(double *ptr) {
  double sum=0.0,ssum;

  for(int i=0;i<ncp;i++) {
    sum+=ptr[clist[me][i]];
  }

  MPI_Allreduce(&sum,&ssum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  return ssum;
}

double DriftSolver::GetAverage(double *ptr) {
  return GetSum(ptr)/occup;
}

void DriftSolver::Dump(std::ostream &ofl) {
  if(me==0) {
    ofl<<"# Drift Solver map file\n"
       <<"# (c) 2012, Yudi Rosandi\n"
       <<"# rosandi@gmail.com\n\n";

    ofl<<"dim "<<nx<<" "<<ny<<" "<<nz<<"\n"
       <<"timestep "<<timestep<<"\n"
       <<"spacestep "<<spacestep<<"\n"
       <<"symbol\n";

    for(int i=0;i<nbound;i++) {
      ofl<<bound[i].sym<<" ";
      char l=bound[i].l==FDCELL_PERIODIC?'p':(bound[i].l==FDCELL_OCCUPIED?'+':'-');
      char r=bound[i].r==FDCELL_PERIODIC?'p':(bound[i].r==FDCELL_OCCUPIED?'+':'-');
      char t=bound[i].t==FDCELL_PERIODIC?'p':(bound[i].t==FDCELL_OCCUPIED?'+':'-');
      char b=bound[i].b==FDCELL_PERIODIC?'p':(bound[i].b==FDCELL_OCCUPIED?'+':'-');
      if(nz!=1) {
        char n=bound[i].n==FDCELL_PERIODIC?'p':(bound[i].n==FDCELL_OCCUPIED?'+':'-');
        char f=bound[i].f==FDCELL_PERIODIC?'p':(bound[i].f==FDCELL_OCCUPIED?'+':'-');
        ofl<<l<<r<<t<<b<<n<<f<<"\n";
      } else  ofl<<l<<r<<t<<b<<"\n";

      /*
      ofl<<" # ("<<bound[i].l<<")"<<"("<<bound[i].r<<")"
      <<"("<<bound[i].t<<")"<<"("<<bound[i].b<<")"
      <<"("<<bound[i].n<<")"<<"("<<bound[i].f<<")";
      */

    }

    ofl<<"\n";

    int m=0;

    ofl<<"map\n";
    for(int k=0;k<nz;k++) {
      for(int j=0;j<ny;j++) {
        for(int i=0;i<nx;i++) {
          ofl<<map_symbol(map[m++]);
        }
      }
    }

    ofl<<"\n";
    ofl.flush();
  }

  MPI_Barrier(MPI_COMM_WORLD);

}

void DriftSolver::Dump(const string fname) {
  ofstream fl(fname.c_str());
  Dump(fl);
  fl.close();
}

void DriftSolver::CopyMap(char* cmap) {
  for (int i=0;i<dim;i++) {
    cmap[i]=map_symbol(map[i]);
  }
}

/**
 this function assign external link to a solver. This is called by external program (domain, console)
**/

void DriftSolver::Link(string pn, double* p) {
  link_updated=false; // mark it!
  // check if Link already exist
  for(int i=0;i<(int)linknm.size();i++)
    if(pn==linknm[i]) {
      string msg("Link name already exists: ");
      msg.append(pn);
      throw msg.c_str();
    }

  link.push_back(p);
  linknm.push_back(pn);
}

// converting a real coordinate to cell index
int DriftSolver::coordtocell(double x, double y, double z) {
  int i=floor(x/spacestep);
  int j=floor(y/spacestep);
  int k=floor(z/spacestep);
  return index(i,j,k);
}

// get map from density array
// caller must take care that the array lengths are correct
// the length of 'syms' is used as length of 'limits'
// make sure that the symbols are listed in the model file!!!!

void DriftSolver::DensityMap(
                             double* dens, // density array. length=dim
                             const char* syms, // list of symbols, NULL terminated
                             double* limits, // limits of all symbols, sorted descending
                             const char* notouch // leave this character untouch or NULL
                             ) {

  int ns=strlen(syms);
  int nc=strlen(notouch);
  char cmap[dim];

  for(int i=0;i<dim;i++) {
    char ch='.',och=map_symbol(map[i]);

    for(int s=0;s<nc;s++)
      if(och==notouch[s]) {ch=och;break;}

    if(ch=='.'){
      for(int s=0;s<ns;s++)
        if(dens[i]>limits[s]) {ch=syms[s]; break;}
    }

    cmap[i]=ch;
  }

  Remap(NULL,cmap);

}

// convert map to array of double

void DriftSolver::MapToArray(double* dmap) {
  for (int i=0;i<dim;i++) {
    dmap[i]=(double)map[i];
  }
}

void DriftSolver::LinkCoupler(DriftCoupler* c) {
  ncoupler++;
  if(ncoupler>MAX_COUPLER) throw "Number of coupler over limit";
  c->init(this);
  coupler[ncoupler-1]=c;
}

// here calculate an optimum timestep
// by checking all cell against the Courant criteria
// better call this after a Step() call to make sure cap&con are actual

double DriftSolver::CourantAdapt(double factor, bool change_dt) {
  double crit=DBL_MAX;
  double newdt;

  for(int i=0;i<ncp;i++) {
    int idx=clist[me][i];
    if(!map[idx]) continue;
    if(con[idx]==0.0) continue; // this should not happen!!

    double cr=0.5*cap[idx]/con[idx]/h2inv; // dt<=h^2/2.D
    // std::cout<<cap[idx]<<"  "<<con[idx]<<cr<<"\n";
    if(crit>cr)crit=cr;
  }

  // here gather minimum!!

  MPI_Allreduce(&crit,&newdt,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);

  // std::cout<<"crit="<<newdt<<"\n";

  newdt*=factor;
  if(change_dt)timestep=newdt;
  return newdt;
}

void DriftSolver::SetAlgorithm(DriftAlgorithm* al){algo=al;algo->init(this);}
void DriftSolver::ClearAlgorithm(){if(algo) delete algo; algo=NULL;}


void DriftAlgorithm::exec() {

  double *A=solver->A;
  double *T=solver->T;
  double *con=solver->con;
  double *cap=solver->cap;
  double h2inv=solver->h2inv;
  int *map=solver->map;
  int me=solver->me;

  for(int i=0;i<solver->ncp;i++) {

    int idx=solver->clist[me][i];
    if(!map[idx]) continue;

    int l,r,t,b,n,f;
    solver->get_neigcell_index(idx,l,r,t,b,n,f);

    A[idx]=0.0;

    double K;

    if(con[idx]==0.0) continue;

    double crit=0.5*cap[idx]/con[idx]/h2inv; // dt<=h^2/2.D

    if(solver->timestep>crit) {
      static char st[1024];
      sprintf(st,"Courant stability criterion violated: "
              "step=%d time_step=%E crit=%E index=%d cap=%E con=%E h^2=%E",
              solver->step, solver->timestep, crit, idx, cap[idx], con[idx], 1./h2inv);
      // FIXME! strict_courant may belong to this class
      if(solver->strict_courant) throw st;
    }

    if(l>=0) if(map[l]) {
      K=con[l]+con[idx];
      if(K!=0.0) {
        K=con[l]*con[idx]/K;
        A[idx] +=K*(T[l]-T[idx]);
      }
    }

    if(r>=0) if(map[r]) {
      K=con[r]+con[idx];
      if(K!=0.0) {
        K=con[r]*con[idx]/K;
        A[idx] +=K*(T[r]-T[idx]);
      }
    }

    if(t>=0) if(map[t]) {
      K=con[t]+con[idx];
      if(K!=0.0) {
        K=con[t]*con[idx]/K;
        A[idx] +=K*(T[t]-T[idx]);
      }
    }

    if(b>=0) if(map[b]) {
      K=con[b]+con[idx];
      if(K!=0.0) {
        K=con[b]*con[idx]/K;
        A[idx] +=K*(T[b]-T[idx]);
      }
    }

    if(solver->nz!=1) {
      if(n>=0) if(map[n]) {
        K=con[n]+con[idx];
        if(K!=0.0) {
          K=con[n]*con[idx]/K;
          A[idx] +=K*(T[n]-T[idx]);
        }
      }
      if(f>=0) if(map[f]) {
        K=con[f]+con[idx];
        if(K!=0.0) {
          K=con[f]*con[idx]/K;
          A[idx] +=K*(T[f]-T[idx]);
        }
      }
    }

  } /*for i:ncp*/
}
