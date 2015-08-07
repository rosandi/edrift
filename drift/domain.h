/* drift class */

#ifndef __DRIFT_DOMAIN__
#define __DRIFT_DOMAIN__

#include <string>
#include <sstream>
#include <fstream>

using std::string;
using std::ifstream;
using std::istringstream;

namespace edrift {
  
  class DriftDomain {
  public:
    string name;
    class DriftSolver* solver;

#ifdef DRIFT_CONSOLE
    std::ofstream output;
#endif
    
    string src;
    string cap;
    string con;
    
    DriftDomain(string infile, string rename="");
    ~DriftDomain();

    void RestoreCon();
    void RestoreCap();
    void RestoreSource();
    
  protected:
    void create_from_file(const char*);

  };
  
}

#endif

