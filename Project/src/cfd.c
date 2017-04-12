#include "cfd.h"


struct _problem* init_problem(){
  return malloc(sizeof(struct _problem));
}

struct _problem* init_problem_physical(struct _problem* Problem, double CFL, double L, double H, double Ls, double Hs, double h, double dt){
  (*Problem).CFL = CFL;
  (*Problem).L   = L;
  (*Problem).H   = H;
  (*Problem).Ls  = Ls;
  (*Problem).Hs  = Hs;
  (*Problem).h   = h;
  (*Problem).dt  = dt;
  return Problem;
}

struct _problem* init_problem_numerical(struct _problem* Problem){
  (*Problem).Nx = (*Problem).H/(*Problem).h;
  (*Problem).Ny = (*Problem).L/(*Problem).h;
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

  free(Problem);
}
