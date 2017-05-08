/*
 *  CFD - Linux Version of the Main
 *  by Son Tran
 */

#include <stdio.h>
#include <stddef.h>
#include <time.h>
//#include <sys/sysinfo.h>
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
  fprintf(stderr, "Problem occured at t= %d\n",(int) ( (*Problem).t * 100 ));
  fprintf(stderr, "Saving current data under  _%d_\n",(int) ( (*Problem).t * 100 ));
  print_problem_data(Problem);
  free_problem_vector_domain(Problem);
  exit(-1);
}

void print_problem_data(struct _problem* Problem){
  char buff_omega_name[50];
  char buff_psi_name[50];
  char buff_u_name[50];
  char buff_v_name[50];
  sprintf(buff_omega_name ,"data/CFD_omega_%d.txt",(int) ( (*Problem).t * 100 ));
  sprintf(buff_psi_name   ,"data/CFD_psi_%d.txt"  ,(int) ( (*Problem).t * 100 ));
  sprintf(buff_u_name     ,"data/CFD_u_%d.txt"    ,(int) ( (*Problem).t * 100 ));
  sprintf(buff_v_name     ,"data/CFD_v_%d.txt"    ,(int) ( (*Problem).t * 100 ));
  FILE* file_omega = fopen(buff_omega_name,"w");
  FILE* file_psi   = fopen(buff_psi_name  ,"w");
  FILE* file_u     = fopen(buff_u_name    ,"w");
  FILE* file_v     = fopen(buff_v_name    ,"w");

  if(file_omega == NULL || file_psi == NULL || file_u == NULL || file_v == NULL){ fprintf(stderr,"File error\n"); exit(1);}
    for(int j = 0 ; j < (*Problem).Ny ; j++){
      for(int i = 0; i < (*Problem).Nx ; i++){
        fprintf(file_omega,"%5.16f "   ,(*Problem).omega[i][j]);
        fprintf(file_psi  ,"%5.16f "   ,(*Problem).psi  [i][j]);
        fprintf(file_u    ,"%5.16f "   ,(*Problem).u    [i][j]);
        fprintf(file_v    ,"%5.16f "   ,(*Problem).v    [i][j]);
    }
    fprintf(file_omega ,"\n");
    fprintf(file_psi   ,"\n");
    fprintf(file_u     ,"\n");
    fprintf(file_v     ,"\n");
  }
  fclose(file_omega);
  fclose(file_psi);
  fclose(file_u);
  fclose(file_v);
}

int main(int argv,char* argc[]){
  // Flow specification
  double nu    =   1e-6;

  // Numerical parameter
  double CFL   =   1.0 ;
  double tau   =   0.1 ;
  double Rey   =   10  ;
  double h     =   0.01;

  // Domain parameter
  double L     =   4.0 ;
  double H     =   1.0 ;
  double Q0    =   Rey*nu;
  double Um    =   Q0/H;
  double Ls    = L/4.0 ;
  double Hs    = H/2.0 ;
  // Computation parameter
  //double Rey_h =
  //double r     =   nu*dt/(h*h);

  double dt    =  0.1;// CFL*h/Um;
  double tmax  =  2.0;// tau*L/Um;

  double phi   =   1.98;
  double tol   =   0.00001;




  fprintf(stderr, "CFD - Oscillating flow in a channel with abrupt step\n");
  fprintf(stderr, "by S. Tran - J. Demey \n");
  fprintf(stderr, "CFL: %f Tau: %f Reynold: %f\n",CFL,tau,Rey);
  fprintf(stderr, "Average Velocity: %f Grid size h: %f Timestep: %f\n",Um,h,dt);
  struct _problem* Problem = init_problem();
  init_problem_physical(Problem,CFL,L,H,Ls,Hs,h,dt,tmax,Q0,tol,nu);
  init_problem_numerical(Problem,phi);
  init_problem_map(Problem);
  init_problem_vector_domain(Problem);
  init_problem_poiseuille(Problem);
  boundary_omega_update(Problem);
  boundary_psi_update(Problem,functionQ);
  //poisson_inner_psi_iterator(Problem);

  print_problem_data(Problem);
  //boundary_u_v_in_out_set(Problem);
  //print_problem_data(Problem);
  // ---Code Benchmarking-------
  struct timespec start, finish;
  double elapsed;
  clock_gettime(CLOCK_MONOTONIC, &start);
  // ---------------------------
  //poisson_inner_psi_iterator(Problem);
  //first_time_integration(Problem);
  integration_omega(Problem);
  // ---Code Benchmarking-------
  clock_gettime(CLOCK_MONOTONIC, &finish);
  elapsed = (finish.tv_sec - start.tv_sec);
  elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
  // ---------------------------
  //test_omega_domainFill (Problem);

  //test_psi_boundaryFill(Problem,test_Qfunc_const);
  print_problem_data(Problem);
  printf("Done \n");


  printf("Time Elapsed: %f s\n",elapsed);
  free_problem_vector_domain(Problem);
}