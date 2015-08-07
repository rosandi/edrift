/* drift class */

#ifndef __DRIFT_BUFFER__
#define __DRIFT_BUFFER__

#include <string>
#include <vector>
#include <fstream>

using std::string;
using std::vector;
using std::ofstream;

namespace edrift {

  class Buffer {
    bool ref; // only refference
    double *prox,*proy,*proz;
    // extra field for dump. Now osed only for vtr
    vector<string> extname;
    vector<Buffer*> extfield;
    
  public:
    string  name;
    bool empty; // this var name is ambigue
    double* ptr;
    int nx,ny,nz;
    
    Buffer();
    Buffer(string vname, int dx, int dy, int dz);
    void alloc(int dx, int dy, int dz);
    void refto(double*,int,int,int);
    void refto(Buffer* bref);
    void refto(class DriftDomain* cref, string vec="cell");
    ~Buffer();
    int size();
    double& value(int ix, int iy, int iz);
    void dump_nz(std::ostream& ofl, int slab);
    void dump_ny(std::ostream& ofl, int slab);
    void dump_nx(std::ostream& ofl, int slab);
    void dump(std::ostream& ofl, int slab, const char norm);
    void dump(std::ostream& ofl, const char norm);
    void dump(std::ostream& ofl, const char norm, int ta, int tb);
    void dump(std::ostream& ofl);
    void dumpxyz(std::ostream& ofl); // xyz sparse
    void dumpcsv(std::ostream& ofl); // comma separated
    void dumpvtr(std::ostream& ofl); // vtk rectilinear
    void dumpvtr(string fname);
    
    void set_extrafield(string name, Buffer& B);
    void clear_extrafield();
    double* prof_z();
    double* prof_y();
    double* prof_x();
    
    /**
     Calculate lateral average (a profile) in axis direction.
     caller is responsible to free (using delete[]) the returned pointer
     */    
    
    double* prof(const char axis) {
      if(axis=='x') return prof_x();
      if(axis=='y') return prof_y();
      if(axis=='z') return prof_z();
      return NULL;
    }

    Buffer& add(Buffer& A);
    Buffer& add(double v);
    Buffer& sub(Buffer& A);  
    Buffer& sub(double v);
    Buffer& mul(Buffer& A);
    Buffer& mul(double v);
    Buffer& div(Buffer& A);
    Buffer& div(double v);
    Buffer& copy(Buffer& A);
    Buffer& fill(double v);
    Buffer& fill(class omd::TableReader& tab, bool derive);
    Buffer& fill(string sym, double delta=1.0, int slab=0, int nslab=-1); 
  };
  
}

#endif
