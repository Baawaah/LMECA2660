/*
 *  CFD - Linux Version of the Main
 *  by Son Tran
 */

#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <sys/sysinfo.h>
#include <unistd.h>
#include "cfd.h"


void deadstop_exit(struct _problem* Problem){
  fprintf(stderr, "[DEADSTOP EXIT] \n");
  fprintf(stderr, "case diagnostic DEADSTOP  \n");
  fprintf(stderr, "    Binary Value for each problem\n");
  fprintf(stderr, "       Mesh Reynold u,v   - 1 \n");
  fprintf(stderr, "       Mesh Reynold omega - 2 \n");
  fprintf(stderr, "       Beta - CFL         - 4 \n");
  fprintf(stderr, "Emergency Exit Procedure Initiate\n");
  fprintf(stderr, "Problem occured at tau = %d\n",(int) ( (*Problem).tau * 100 ));
  fprintf(stderr, "Saving current data under  _%d_\n",(int) ( (*Problem).tau * 100 ));
  print_problem_data(Problem);
  free_problem_vector_domain(Problem);
  exit(-1);
}

void print_problem_data(struct _problem* Problem){
  char buff_omega_name[50];
  char buff_psi_name[50];
  char buff_u_name[50];
  char buff_v_name[50];
  char buff_R_name[50];
  sprintf(buff_omega_name ,"data/CFD_omega_%d.txt",(int) ( (*Problem).tau * 100 ));
  sprintf(buff_psi_name   ,"data/CFD_psi_%d.txt"  ,(int) ( (*Problem).tau * 100 ));
  sprintf(buff_u_name     ,"data/CFD_u_%d.txt"    ,(int) ( (*Problem).tau * 100 ));
  sprintf(buff_v_name     ,"data/CFD_v_%d.txt"    ,(int) ( (*Problem).tau * 100 ));
  sprintf(buff_R_name     ,"data/CFD_R_%d.txt"    ,(int) ( (*Problem).tau * 100 ));
  FILE* file_omega = fopen(buff_omega_name,"w");
  FILE* file_psi   = fopen(buff_psi_name  ,"w");
  FILE* file_u     = fopen(buff_u_name    ,"w");
  FILE* file_v     = fopen(buff_v_name    ,"w");
  FILE* file_R     = fopen(buff_R_name    ,"w");

  if(file_omega == NULL || file_psi == NULL || file_u == NULL || file_v == NULL){ fprintf(stderr,"File error\n"); exit(1);}
    for(int j = 0 ; j < (*Problem).Ny ; j++){
      for(int i = 0; i < (*Problem).Nx ; i++){
        fprintf(file_omega,"%5.16f "   ,(*Problem).omega[i][j]);
        fprintf(file_psi  ,"%5.16f "   ,(*Problem).psi  [i][j]);
        fprintf(file_u    ,"%5.16f "   ,(*Problem).u    [i][j]);
        fprintf(file_v    ,"%5.16f "   ,(*Problem).v    [i][j]);
        fprintf(file_R    ,"%5.16f "   ,(*Problem).R_res[i][j]);
    }
    fprintf(file_omega ,"\n");
    fprintf(file_psi   ,"\n");
    fprintf(file_u     ,"\n");
    fprintf(file_v     ,"\n");
    fprintf(file_R     ,"\n");
  }
  fclose(file_omega);
  fclose(file_psi);
  fclose(file_u);
  fclose(file_v);
  fclose(file_R);
}

int main(int argv,char* argc[]){
  // Flow specification
  double nu    =   1e-6;

  // Numerical parameter
  double CFL        =   0.2 ;
  double tau        =   1.0 ;
  double Rey        =   100 ;
  double h          =   0.02;
  double t_snapshot =   0.25;

  // Domain parameter
  double L     =   4.0 ;
  double H     =   1.0 ;
  double Q0    =   Rey*nu;
  double Um    =   Q0/H;
  double Ls    =   L/4.0 ;
  double Hs    =   H/2.0 ;
  // Computation parameter
  //double Rey_h =
  //double r     =   nu*dt/(h*h);

  double dt    =  CFL*h/Um;
  double tmax  =  tau*L/Um;

  double phi   =   1.98;
  double tol   =   1e-4;




  fprintf(stderr, "CFD - Oscillating flow in a channel with abrupt step\n");
  fprintf(stderr, "by S. Tran - J. Demey \n");
  fprintf(stderr, "CFL: %f Tau: %f Reynold: %f\n",CFL,tau,Rey);
  fprintf(stderr, "Average Velocity: %f Grid size h: %f Timestep: %f\n",Um,h,dt);

  struct _problem* Problem = init_problem();
  init_problem_physical(Problem,CFL,L,H,Ls,Hs,h,dt,tmax,Q0,tol,nu,Um);
  init_problem_numerical(Problem,phi,t_snapshot);
  init_problem_map(Problem);
  init_problem_vector_domain(Problem);
  init_problem_poiseuille(Problem);
  boundary_psi_update(Problem,functionQ);
  poisson_inner_psi_iterator(Problem);
  boundary_omega_update(Problem);
  //boundary_omega_dwdx_update(Problem);
  inner_u_v_compute(Problem);

  fprintf(stderr, "NX: %d NY: %d NLs: %d NHs: %d \n",(*Problem).Nx,(*Problem).Ny,(*Problem).NLs,(*Problem).NHs);
  print_problem_data(Problem);
  printf("Simulation Starting\n");
  // ---Code Benchmarking-------
  struct timespec start, finish;
  double elapsed;
  clock_gettime(CLOCK_MONOTONIC, &start);
  // ---------------------------

  //integration_omega(Problem);

  // ---Code Benchmarking-------
  clock_gettime(CLOCK_MONOTONIC, &finish);
  elapsed = (finish.tv_sec - start.tv_sec);
  elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
  // ---------------------------

  //print_problem_data(Problem);
  printf("Simulation Done\n");


  printf("Time Elapsed: %f s\n",elapsed);
  printf("Final File save under _%d_\n",(int) ( (*Problem).tau * 100 ));
  free_problem_vector_domain(Problem);
}
