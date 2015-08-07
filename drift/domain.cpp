/* drift class */

#include "domain.h"
#include "drift.h"
#include "omd/omdtool.h"

using namespace omd;
using namespace edrift;

string getall(istringstream& ss) {
  string str;
  while(ss.good()) {
    string st;
    if(ss>>st) str.append(st+" "); // space separated
  }
  return str;
}

string read_block(ifstream &mfile) {
  char sline[4096];
  string sblock;
  string tok;

  // a block ends with empty line
  while(mfile.good()) {
    mfile.getline(sline,4096);
    istringstream ss(sline);
    tok="";
    if(!(ss>>tok)) break;
    if(tok[0]=='#') continue;
    sblock.append(sline);
    sblock.append(" "); // append space
  }

  return sblock;
}

// rad map data
string read_map(ifstream &mfile) {
  char sline[4096];
  string sblock;
  string stmp;
  string tok;
  bool inbra=false;
  int mul, nsl;

  while(mfile.good()) {
    mfile.getline(sline,4096);
    istringstream ss(sline);
    if(!(ss>>tok)) break;

    if(tok=="{") {
      if(inbra) throw "no nesting bracket allowed";
      if(!(ss>>nsl)) nsl=1;
      inbra=true;

      sblock.append(stmp);
      stmp.clear();
      continue;

    } else if(tok=="}") {
      for(int i=0;i<nsl;i++) sblock.append(stmp);
      stmp.clear();
      inbra=false;
      continue;
    }

    if(!(ss>>mul)) mul=1;
    for(int i=0;i<mul;i++)
      stmp.append(tok);

  }
  if(inbra) throw "missing closing bracket";
  sblock.append(stmp);
  return sblock;
}


DriftDomain::DriftDomain(string infile, string rename) {
  create_from_file(infile.c_str());
  if(rename.size()) name=rename;
}

DriftDomain::~DriftDomain() {
  delete solver;
}

void DriftDomain::create_from_file(const char* mapfile) {
  ifstream mfile(mapfile);

  if(mfile.fail()) {
    string msg("failed to open file: ");
    msg.append(mapfile);
    throw msg.c_str();
  }

  string map="", symbol="", con_string="", cap_string="", src_string="";
  string cell_string="", var_string="", equ_string="";
  string tab_string="", string_pbc="",coup_string="", link_string="";

  int  nx,ny,nz;
  double dt=1.0,h=1.0;
  char sline[4096];
  bool strict_courant=true;

  while(mfile.good()) {
    mfile.getline(sline,4096);
    istringstream ss(sline);
    string tok;
    ss>>tok;

    // skip empty spaces and comments
    if(tok=="") continue;
    if(tok.at(0)=='#') continue;


    if(tok=="dim") {
      ss >> nx >> ny;
      if(!(ss>>nz)) nz=1; // one dimension
    }
    else if(tok=="timestep") ss>>dt;
    else if(tok=="spacestep") ss>>h;
    else if(tok=="name") ss>>name;
    else if(tok=="periodic") ss>>string_pbc;

    else if(tok=="symbol") {
      symbol=read_block(mfile);
    }

    else if(tok=="map") {
      /*
      string mapfl("");
      if(ss>>mapfl) {
        if(mapfl.at(0)=='#') mapfl="";
      } else throw (string("can not open file: ")+mapfl).c_str();

      if(mapfl=="") map=read_map(mfile);
      else {
        ifstream mfl(mapfl.c_str());
        map=read_map(mfl);
        mfl.close();
      }
       */

      map=read_map(mfile);
    }

    else if(tok=="conductivity") {
      con_string=read_block(mfile);
    }

    else if(tok=="capacity") {
      cap_string=read_block(mfile);
    }

    else if(tok=="source") {
      src_string=read_block(mfile);
    }

    else if(tok=="cell") {
      cell_string=read_block(mfile);
    }

    else if(tok=="coupling") {
      coup_string=read_block(mfile);
    }

    else if(tok=="table") {
      string ftab,scl,c;
      if(!(ss>>ftab)) throw "table name required: [table@]filename";
      tab_string.append(ftab+" ");
      if(ss>>scl>>c) tab_string.append(scl+" "+c+" ");
      else tab_string.append("1.0 1.0 ");
    }

    else if(tok=="set") {
      string expr(remove_char(getall(ss),' '));
      if(expr.find('=')==expr.npos) throw "invalid set command";
      var_string.append(expr+" ");
    }

    else if(tok=="extern") {
      string nm;
      if(ss>>nm) link_string.append(nm+" ");
    }

    else if(tok=="expr") {
      string expr(remove_char(getall(ss),' '));
      if(expr.find('=')==expr.npos) throw "invalid expr command";
      equ_string.append(expr+" ");
    }

    else if(tok=="strict_courant") {
      string sbool;
      if(ss>>sbool) { // no value means true
        if(sbool=="yes") strict_courant=true;
        if(sbool=="true") strict_courant=true;
        if(sbool=="no") strict_courant=false;
        if(sbool=="false") strict_courant=false;
      }
    }

    else {
      string err("unrecognized token: ");
      err.append(tok);
      throw err.c_str();
    }

  }
  mfile.close();

  if(nx<=0 || ny<=0 || nz<=0) throw "no dimension definition in file (dim)";
  if(symbol=="" || map=="") throw "no symbol and/or map defined in file";

  solver=new DriftSolver(nx,ny,nz,symbol.c_str(),map.c_str());
  solver->SetAlgorithm(new DriftAlgorithm);

  solver->strict_courant=strict_courant;
  solver->SetTimestep(dt);
  solver->SetSpacestep(h);

  int pbc=0;
  if(string_pbc.find('x')!=string_pbc.npos) pbc|=1;
  if(string_pbc.find('y')!=string_pbc.npos) pbc|=2;
  if(string_pbc.find('z')!=string_pbc.npos) pbc|=4;

  if(pbc) solver->ForcePeriodic(pbc);

  src=src_string;
  con=con_string;
  cap=cap_string;

  // this function checks everything before running simulation
  solver->CheckStrings();

}
