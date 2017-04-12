/*
 *  CFD - Linux Version of the Main
 *  by Son Tran
 */

#include <stdio.h>
#include <time.h>
#include <sys/sysinfo.h>
#include "cfd.h"

int main(int argv,char* argc[]){
  double CFL =   1.0;
  double L   =   1.0;
  double H   =   0.5;
  double h   =   0.1;
  double dt  =   0.1;
  double Ls  = L/2.0;
  double Hs  = H/2.0;

  struct _problem* Problem = init_problem_vector_domain(init_problem_numerical(init_problem_physical(init_problem(),CFL,L,H,Ls,Hs,h,dt)));
  free_problem_vector_domain(Problem);
 printf("Simple Check");
}
