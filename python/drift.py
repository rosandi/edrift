import sys,traceback,types
from ctypes import *

class drift:
    def __init__(self,nx=0,ny=0,nz=0,driftsym=None,driftmap=None,mapfile=None):
        try:
            self.lib = CDLL("libdrift.so",RTLD_GLOBAL)
        except:
            ertype,value,tb = sys.exc_info()
            traceback.print_exception(ertype,value,tb)
            raise OSError,"Could not load EDRIFT dynamic library"

        self.drift=c_void_p()
        if(mapfile):
            self.read_mapfile(mapfile)
        else:
            self.dim=nx*ny*nz;
            self.lib.drift_create(nx,ny,nz,driftsym,driftmap,byref(self.drift))
            self.created=1

    def read_mapfile(self,mapfile):
        print 'reading mapfile ',mapfile
        nx=0;ny=0;nz=0
        strsym=""
        strmap=""
        try:
            fl=open(mapfile,'r')
            for line in fl:
                line=line.rstrip().split(' ')
                if (line[0]=='dim'):
                    nx=int(line[1])
                    ny=int(line[2])
                    nz=int(line[3])
                    continue

                if (line[0]=='symbol'):

                    for ss in fl:
                        if ss=='\n': break
                        strsym+=ss.rstrip()+' '
                    continue

                if (line[0]=='map'):
                    for ss in fl:
                        if ss=='\n': break
                        strmap+=ss.rstrip().replace(' ','')
                    continue

                # other tokens are ignored

            fl.close()

            if (nx==0 or ny==0 or nz==0):
                raise OSError,"MAPFILE read failure"

            self.dim=nx*ny*nz
            self.lib.drift_create(nx,ny,nz,strsym,strmap,byref(self.drift))
            self.created=1

        except:
            ertype,value,tb = sys.exc_info()
            traceback.print_exception(ertype,value,tb)
            raise OSError,"MAPFILE read failure"

    def __del__(self):
        if self.drift and self.created: self.lib.drift_destroy(self.drift)

    def timestep(self,dt):
        self.lib.drift_set_timestep(self.drift,c_double(dt))

    def spacestep(self,dx):
        self.lib.drift_set_spacestep(self.drift,c_double(dx))

    def dimension(self):
        nx=self.lib.drift_getdim(self.drift,0)
        ny=self.lib.drift_getdim(self.drift,1)
        nz=self.lib.drift_getdim(self.drift,2)
        return [nx,ny,nz]

    def algorithm(self,algoname):
        self.lib.drift_set_algorithm(self.drift,algoname)

    def step(self,n=1):
        self.lib.drift_step(self.drift,n)

    def cell(self):
        self.lib.drift_cell_vec.restype = POINTER(c_double)
        ptr=self.lib.drift_cell_vec(self.drift);
        return ptr

    def source(self):
        self.lib.drift_source_vec.restype = POINTER(c_double)
        ptr=self.lib.drift_source_vec(self.drift)
        return ptr

    def ext(self):
        self.lib.drift_ext_vec.restype = POINTER(c_double)
        ptr=self.lib.drift_ext_vec(self.drift);
        return ptr

    def con(self):
        self.lib.drift_con_vec.restype = POINTER(c_double)
        ptr=self.lib.drift_con_vec(self.drift);
        return ptr

    def cap(self):
        self.lib.drift_cap_vec.restype = POINTER(c_double)
        ptr=self.lib.drift_cap_vec(self.drift);
        return ptr

    def setcell(self, expr, cells="*"):
        self.lib.drift_setcell(self.drift,c_double(expr),cells)

    def setsource(self, expr, cells="*"):
        self.lib.drift_setsource(self.drift,c_double(expr),cells)

    def setext(self, expr, cells="*"):
        self.lib.drift_setext(self.drift,c_double(expr),cells)

    def setcap(self,expr, cells="*"):
        self.lib.drift_setcap(self.drift,c_double(expr),cells)

    def setcon(self, expr, cells="*"):
        self.lib.drift_setcon(self.drift,c_double(expr),cells)

    def dump(self, fname):
        self.lib.drift_dump(self.drift,fname)

    def collect(self,ptr):
        self.lib.drift_collect(self.drift,ptr)
        i=0
        a=[]
        while (i<self.dim):
            a.append(ptr[i])
            i+=1
        return a

    def link_equhead(self,ptr,equ='0',symbol='*',trange=[0,sys.maxint]):
        self.lib.drift_equ_head(self.drift,ptr,equ,symbol,c_double(trange[0]),c_double(trange[1]))

    def link_equtail(self,ptr,equ='0',symbol='*',trange=[0,sys.maxint]):
        self.lib.drift_equ_tail(self.drift,ptr,equ,symbol,c_double(trange[0]),c_double(trange[1]))

    def link_source(self,e0,table=None,pulsewidth=0):
        if table:
            self.lib.drift_source_table(self.drift,e0,table)
        else:
            self.lib.drift_source_pulse(self.drift,e0,pulsewidth)

    def link_prolb(self,e0,lamb,n0,nmap,pulsewidth=0,tname=None):
        if tname:
            self.lib.drift_prolb_table(self.drift,e0,lamb,n0,nmap,tname)
        else:
            self.lib.drift_prolb_table(self.drift,e0,lamb,n0,nmap,pulsewidth)

    def link_ttmcoupling(self,te,ta,coupconst=0,couptable=None):
        if couptable:
            self.lib.drift_ttm_coupling_table(self.drift, couptable,te,ta)
        else:
            self.lib.drift_ttm_coupling(self.drift,coupconst,te,ta)

    def link_capteta(self,k0,te,ta):
        self.lib.drift_teta(self.drift,k0,te,ta)

    def link_tablefill(self,tname,xarray,target):
        self.lib.drift_table(self.drift,tname,xarray,target)
