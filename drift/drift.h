
#ifndef __DRIFT_H__
#define __DRIFT_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <omd/treader.h>
#include <unistd.h>

// For 2D problems...

using std::string;
using std::vector;
using std::ifstream;
using std::istringstream;

namespace edrift {

#define FDCELL_EMPTY 0
#define FDCELL_OCCUPIED 1
#define FDCELL_PERIODIC 2
#define MAX_COUPLER 10
#define FD_MAX_VAR 50

  struct CellBoundary {
    char sym;
    // left right top bottom near far
    int  l,r,t,b,n,f;
  };

  class DriftSolver {

    class DriftAlgorithm* algo;
    bool link_updated; // unset by Link(), set by UpdateLink()

  public:

    double*  A; // cell value array
    double*  T; // temporary cell array
    double*  S; // source array
    double*  Ext; // external array (coupling, etc)

    double* dE;
    double* dS;
    double* dExt;
    double* sumE;
    double* sumS;
    double*  con;
    double*  cap;

    int ncoupler;
    class DriftCoupler* coupler[MAX_COUPLER];

    vector<double*> link;
    vector<string> linknm;

    int me; // mpi rank
    int nproc; // number of proc
    int occup; // total # of occupied cells
    int ncp;
    int ncell;
    int *lim; // hi/lo limits of all procs
    int* counts;
    int* displ;
    int* sendnum;
    int* recvnum;
    int** nlist; // index to send
    int** glist; // index to receive
    double** ndata; // data to send
    double** gdata; // data to receive
    int** clist; // occupied cell list (not empty) of all procs

    int dim,nx,ny,nz;
    int forced_periodic;
    int step;
    double timestep;
    double spacestep;
    double h2inv; // = 1/spacestep^2
    bool strict_courant;

    // map contains index to the coeficients
    int* map;

    int nbound;
    CellBoundary* bound;

    int search_symbol(char c);
    char map_symbol(int num);

    void get_cellcoord(int idx, int& i, int& j, int& k){
      int xy=nx*ny;
      int a=idx%xy;
      k=idx/xy;
      j=a/nx;
      i=a%nx;
    }

    void initiate_mem();
    void initiate_comm();
    void sync_data();
    void get_neigcell_index(int,int&,int&,int&,int&,int&,int&);
    void send_receive_list();
    void select_ncell(int,int,char**);
    bool mycell(int,int,int);

  public:
    int index(int i, int j, int k) {
      int idx=((k*nx*ny)+(i+j*nx));
      if(idx<0||idx>=dim) return -1;
      return idx;
    }

    void equation_setup();

    vector<string> tabname;

    string name;

    DriftSolver(int, int, int, string, string);
    virtual ~DriftSolver();

    void communicate(double* ptr);
    void collect(double* ptr);
    void collect(); // collect all

    void GetIndexList(const char*, vector<int>&);
    void AssignArray(double*,int,int,int,double); // use grid indices
    void AssignArray(double*,const char*,double); // use symbol(s)

    double* GetSourceVector() {return S;}
    double* GetExtVector() {return Ext;}
    double* GetPostSourceVector() {return dS;}
    double* GetPostExtVector() {return dExt;}
    double* GetVector(){return A;}
    double* GetConductivity(){return con;}
    double* GetCapacity(){return cap;}
    double* GetDelta(){return dE;}
    double* GetSourceTally(){return sumS;}
    double* GetExtTally(){return sumE;}
    double GetTimestep(){return timestep;}
    double GetSpacestep(){return spacestep;}
    int  GetStep(){return step;}
    int GetDim(){return dim;}
    void GetDim(int& _nx, int& _ny, int& _nz){_nx=nx;_ny=ny;_nz=nz;}
    char GetMapSymbol(int idx){return map_symbol(map[idx]);}

    double GetSum(double*);
    double GetAverage(double*);

    void LinkCoupler(class DriftCoupler*);
    void UnlinkCoupler(){ncoupler=0;} // this unlink all coupler! Relink any needed coupler.

    void CopyMap(char*); // convert map array to character map

    void SetSymbolMap(string str);

    void SetMap(const char*);
    void Remap(const char*,const char*);   // change completely cell map

    void EnableDelta();
    void EnableTally();
    void KeepSource();

    void SetCell(const char *st, double val){AssignArray(A,st,val);}
    void SetExt(const char *st, double val){AssignArray(Ext,st,val);}
    void SetSource(const char *st, double val){AssignArray(S,st,val);}
    void SetConductivity(const char *st, double val){AssignArray(con,st,val);}
    void SetCapacity(const char *st, double val){AssignArray(cap,st,val);}

    void SetTimestep(double dt);
    void SetSpacestep(double h);

    void Step(int n);
    void Dump(std::ostream& ofl);
    void Dump(const string fname);
    void ForcePeriodic(int bitmask) {forced_periodic=bitmask;} // (00000ZYX)b

    int GetNCP(){return ncp;} // # of local cell
    int* GetCellList(){return clist[me];}

    void Link(string, double*);

    int coordtocell(double x, double y, double z);

    void DensityMap(
                    double* dens,
                    const char* syms,
                    double* limits,
                    const char* notouch=NULL);

    void MapToArray(double* dmap);

    void update();

    // FIXME! to check all cell and give new timestep
    // factor: the factor of the criterion
    // change_dt: only calculate of put new value on timestep
    double CourantAdapt(double factor=1.0, bool change_dt=false);

    void SetAlgorithm(DriftAlgorithm* al);
    void ClearAlgorithm();

  };

  class DriftAlgorithm {
  public:
    DriftSolver* solver;
    DriftAlgorithm(){solver=NULL;}
    virtual ~DriftAlgorithm(){}
    virtual void init(DriftSolver* slv){solver=slv;}
    virtual void exec();
  };

  // to be inherited by any coupling algorithms

  class DriftCoupler {
  public:
    DriftSolver* solver;
    DriftCoupler(){solver=NULL;}
    virtual ~DriftCoupler(){}
    virtual void init(DriftSolver* slv){solver=slv;}
    virtual void update(){}
    virtual void head(){}
    virtual void tail(){}
  };

}

#endif
