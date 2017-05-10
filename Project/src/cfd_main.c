#include "cfd.h"



struct _problem* init_problem(){
  return malloc(sizeof(struct _problem));
}

void init_problem_physical(struct _problem* Problem, double CFL, double L, double H, double Ls, double Hs, double h, double dt, double tmax, double Q0, double tol, double nu){
  (*Problem).CFL   = CFL;
  (*Problem).L     = L;
  (*Problem).H     = H;
  (*Problem).Ls    = Ls;
  (*Problem).Hs    = Hs;
  (*Problem).h     = h;
  (*Problem).dt    = dt;
  (*Problem).t     = 0;
  (*Problem).tmax  = tmax;
  (*Problem).Q0    = Q0;
  (*Problem).tol   = tol;
  (*Problem).nu    = nu;
}

//comment

void init_problem_numerical(struct _problem* Problem, double phi){
  (*Problem).Nx    = (*Problem).L   /(*Problem).h;
  (*Problem).Ny    = (*Problem).H   /(*Problem).h;
  (*Problem).Ntime = (*Problem).tmax/(*Problem).dt;
  (*Problem).phi   = phi;
}

void init_problem_map(struct _problem* Problem){
  (*Problem).imax_map = calloc((*Problem).Nx,sizeof(int));
  for(int i = 0 ; i < (*Problem).Nx ; i++){
      if( (*Problem).h*i < (*Problem).Ls ) (*Problem).imax_map[i] = ((*Problem).H - (*Problem).Hs)/(*Problem).h ;
      else (*Problem).imax_map[i] = (*Problem).Ny;
  }
  (*Problem).NLs = (*Problem).Ls/(*Problem).h;
  (*Problem).NHs = (*Problem).Hs/(*Problem).h;
}

void init_problem_vector_domain(struct _problem* Problem){
  (*Problem).omega = (double**) malloc((*Problem).Nx*sizeof(double*));
  (*Problem).psi   = (double**) malloc((*Problem).Nx*sizeof(double*));
  (*Problem).psi_s = (double**) malloc((*Problem).Nx*sizeof(double*));
  (*Problem).u     = (double**) malloc((*Problem).Nx*sizeof(double*));
  (*Problem).v     = (double**) malloc((*Problem).Nx*sizeof(double*));
  (*Problem).f_old = (double**) malloc((*Problem).Nx*sizeof(double*));
  (*Problem).R_res = (double**) malloc((*Problem).Nx*sizeof(double*));
  for(int i=0 ; i < (*Problem).Nx ; i++){
    (*Problem).omega[i] = (double*) calloc((*Problem).Ny,sizeof(double));
    (*Problem).psi[i]   = (double*) calloc((*Problem).Ny,sizeof(double));
    (*Problem).psi_s[i] = (double*) calloc((*Problem).Ny,sizeof(double));
    (*Problem).u[i]     = (double*) calloc((*Problem).Ny,sizeof(double));
    (*Problem).v[i]     = (double*) calloc((*Problem).Ny,sizeof(double));
    (*Problem).f_old[i] = (double*) calloc((*Problem).Ny,sizeof(double));
    (*Problem).R_res[i] = (double*) calloc((*Problem).Ny,sizeof(double));
  }
}

void init_problem_poiseuille(struct _problem* Problem){
  for(int i = 1; i < (*Problem).Nx-1;i++ ){
    for(int j = 1; j < (*Problem).NHs-1; j++ ){
      double eta =  ( ((*Problem).Hs/2.0) - (j)*(*Problem).h )/((*Problem).Hs/2.0);
      (*Problem).u[i][j]     = scalar_u_v_poiseuille(Problem,eta);
      (*Problem).omega[i][j] = -( (*Problem).u[i][j] - (*Problem).u[i][j-1]) /(*Problem).h ;
    }
  }
}


void free_problem_vector_domain(struct _problem* Problem){
  for(int i=0 ; i < (*Problem).Nx ; i++){
    free((*Problem).omega[i]);
    free((*Problem).psi[i]);
    free((*Problem).psi_s[i]);
    free((*Problem).u[i]);
    free((*Problem).v[i]);
    free((*Problem).f_old[i]);
    free((*Problem).R_res[i]);
  }
  free((*Problem).omega);
  free((*Problem).psi);
  free((*Problem).psi_s);
  free((*Problem).u);
  free((*Problem).v);
  free((*Problem).f_old);
  free((*Problem).R_res);
  free((*Problem).imax_map);

  free(Problem);
}
