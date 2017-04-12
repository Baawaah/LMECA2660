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

  for(int i = 0; i < (*Problem).Nx ; i++){
    for(int j = 0 ; j < (*Problem).Ny ; j++){
        fprintf(file_omega,"%f "   ,(*Problem).omega[i][j]);
        fprintf(file_psi  ,"%f "   ,(*Problem).psi  [i][j]);
        fprintf(file_u    ,"%f "   ,(*Problem).u    [i][j]);
        fprintf(file_v    ,"%f "   ,(*Problem).v    [i][j]);
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
  double CFL  =   1.0;
  double L    =   1.0;
  double H    =   0.5;
  double h    =   0.1;
  double dt   =   0.1;
  double Ls   = L/2.0;
  double Hs   = H/2.0;
  double tmax =   1.0;

  struct _problem* Problem = init_problem_vector_domain(init_problem_map(init_problem_numerical(init_problem_physical(init_problem(),CFL,L,H,Ls,Hs,h,dt,tmax))));
  print_problem_data(Problem);
  free_problem_vector_domain(Problem);
  printf("Check \n");
}
