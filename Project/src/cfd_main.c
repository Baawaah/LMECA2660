#include "cfd.h"



struct _problem* init_problem(){
  return malloc(sizeof(struct _problem));
}

struct _problem* init_problem_physical(struct _problem* Problem, double CFL, double L, double H, double Ls, double Hs, double h, double dt, double tmax){
  (*Problem).CFL = CFL;
  (*Problem).L   = L;
  (*Problem).H   = H;
  (*Problem).Ls  = Ls;
  (*Problem).Hs  = Hs;
  (*Problem).h   = h;
  (*Problem).dt  = dt;
  (*Problem).t  = 0;
  (*Problem).tmax  = tmax;
  return Problem;
}

//comment

struct _problem* init_problem_numerical(struct _problem* Problem){
  (*Problem).Nx    = (*Problem).L   /(*Problem).h;
  (*Problem).Ny    = (*Problem).H   /(*Problem).h;
  (*Problem).Ntime = (*Problem).tmax/(*Problem).dt;
  return Problem;
}

struct _problem* init_problem_map(struct _problem* Problem){
  (*Problem).imax_map = calloc((*Problem).Nx,sizeof(int));
  for(int i = 0 ; i < (*Problem).Nx ; i++){
      if( (*Problem).h*i < (*Problem).Ls ) (*Problem).imax_map[i] = ((*Problem).H - (*Problem).Hs)/(*Problem).h ;
      else (*Problem).imax_map[i] = (*Problem).Ny;
  }
  (*Problem).NLs = (*Problem).Ls/(*Problem).h;
  (*Problem).NHs = (*Problem).Hs/(*Problem).h;
  return Problem;
}

struct _problem* init_problem_vector_domain(struct _problem* Problem){
  (*Problem).omega = (double**) malloc((*Problem).Nx*sizeof(double*));
  (*Problem).psi   = (double**) malloc((*Problem).Nx*sizeof(double*));
  (*Problem).u     = (double**) malloc((*Problem).Nx*sizeof(double*));
  (*Problem).v     = (double**) malloc((*Problem).Nx*sizeof(double*));
  for(int i=0 ; i < (*Problem).Nx ; i++){
    (*Problem).omega[i] = (double*) calloc((*Problem).Ny,sizeof(double));
    (*Problem).psi[i]   = (double*) calloc((*Problem).Ny,sizeof(double));
    (*Problem).u[i]     = (double*) calloc((*Problem).Ny,sizeof(double));
    (*Problem).v[i]     = (double*) calloc((*Problem).Ny,sizeof(double));
  }
  return Problem;
}

void free_problem_vector_domain(struct _problem* Problem){
  for(int i=0 ; i < (*Problem).Nx ; i++){
    free((*Problem).omega[i]);
    free((*Problem).psi[i]);
    free((*Problem).u[i]);
    free((*Problem).v[i]);
  }
  free((*Problem).omega);
  free((*Problem).psi);
  free((*Problem).u);
  free((*Problem).v);

  free((*Problem).imax_map);

  free(Problem);
}
