
#ifdef __cplusplus
extern "C" {
#endif
void drift_create(int nx, int ny, int nz, char* sym, char* map, void** ptr);
void drift_destroy(void* ptr);
void drift_set_timestep(void* ptr, double timestep);
void drift_set_spacestep(void* ptr,double delta);
void drift_set_algorithm(void* ptr,const char* algoname);

void drift_step(void* ptr,int nstep);
void drift_dump(void* ptr,const char* fname);

void* drift_cell_vec(void* ptr);
void* drift_source_vec(void* ptr);
void* drift_ext_vec(void* ptr);
void* drift_con_vec(void* ptr);
void* drift_cap_vec(void* ptr);

void drift_setcell(void* ptr, double val, const char* st);
void drift_setsource(void* ptr, double val, const char* st);
void drift_setext(void* ptr, double val, const char* st);
void drift_setcon(void* ptr, double val, const char* st);
void drift_setcap(void* ptr, double val, const char* st);

void drift_collect(void* ptr, double* arr);
double drift_getdim(void* ptr, int axis);

void drift_teta(void* ptr, double k, double* te, double* ta);
void drift_table(void* ptr, const char* tname, double* x, double* y);

void drift_ttm_coupling_table(void* ptr, const char* tname, double* te, double* ta);
void drift_ttm_coupling(void* ptr, double cval, double* te, double* ta);
void drift_prolb(void* ptr, double e0, double lamb, double n0,  double* nmap, double tau);
void drift_prolb_table(void* ptr, double e0, double lamb, double n0, double* nmap, const char* tname);
void drift_source_pulse(void* ptr, double e0,double pwidth);
void drift_source_table(void* ptr, double e0, const char* tname);
void drift_equ_head(void* ptr, double* data, const char* equ, const char* cells, int sta, int sto);
void drift_equ_tail(void* ptr, double* data, const char* equ, const char* cells, int sta, int sto);

#ifdef __cplusplus
}
#endif
