#include "cfd.h"
#include <unistd.h>
#include <stdio.h>
#include <stddef.h>

double functionQ(struct _problem* Problem){
  return  (*Problem).Q0; //(*Problem).Q0*cos(2*M_PI*t);
}

double scalar_rhs_conv(struct _problem* Problem, int i, int j){
    double**omega_loc = (*Problem).omega;
    double domdx,domdy;

    domdx = (*Problem).u[i][j]*( omega_loc[i+1][j]   -  omega_loc[i-1][j]   )/(2.0*(*Problem).h);
    domdy = (*Problem).v[i][j]*( omega_loc[i]  [j+1] -  omega_loc[i]  [j-1] )/(2.0*(*Problem).h);

    return (domdy + domdx);

}

double scalar_rhs_diff(struct _problem* Problem, int i, int j){
    double**omega_loc = (*Problem).omega;
    double dom2dx2,dom2dy2 ;

    dom2dx2 = (*Problem).nu*( omega_loc[i+1][j]  - 2.0*omega_loc[i][j] + omega_loc[i-1][j]   )/pow((*Problem).h,2.0);
    dom2dy2 = (*Problem).nu*( omega_loc[i]  [j+1]- 2.0*omega_loc[i][j] + omega_loc[i]  [j-1] )/pow((*Problem).h,2.0);

    return (dom2dy2 + dom2dx2);

}

void first_iteration_omega(struct _problem* Problem){
	double dom_old_conv,dom_old_diff;
   (*Problem).t = (*Problem).t + (*Problem).dt;
    for(int i=1; i<(*Problem).Nx-1; i++){
        for(int j = (*Problem).imap[i]; j < (*Problem).Ny-1; j++){
            dom_old_conv = scalar_rhs_conv(Problem,i,j);
            dom_old_diff = scalar_rhs_diff(Problem,i,j);
            (*Problem).omega[i][j]= (*Problem).omega[i][j] - (*Problem).dt*dom_old_conv + (*Problem).dt*dom_old_diff;
            (*Problem).f_old[i][j]= dom_old_conv;
        }
    }

    poisson_inner_psi_iterator(Problem);
    boundary_omega_update(Problem);
    boundary_omega_dwdx_update(Problem);
    inner_u_v_compute(Problem);
}

void integration_omega(struct _problem* Problem){
    int check = 0;
    double dom_conv,dom_diff;
    first_iteration_omega(Problem);
    print_problem_data(Problem);

    for(int k=1; k < (*Problem).Ntime; k++ ){
      (*Problem).t = (*Problem).t + (*Problem).dt;
      // ## Calcule Omega+1
      for(int i=1; i<(*Problem).Nx-1; i++){
          for(int j = (*Problem).imap[i]+1; j < (*Problem).Ny-1; j++){
              dom_conv=scalar_rhs_conv(Problem,i,j);
              dom_diff=scalar_rhs_diff(Problem,i,j);
              (*Problem).omega[i][j] = (*Problem).omega[i][j] - (*Problem).dt*0.5*(3.0*dom_conv - (*Problem).f_old[i][j] )+ (*Problem).dt*dom_diff ;
              (*Problem).f_old[i][j] = dom_conv;

              check = diagnose_check(Problem,i,j);
              //if(check != 0 ){ fprintf(stderr, "[DEADSTOP] Diagnose Check: %d\n",check); deadstop_exit(Problem);}
          }
      }

      boundary_omega_update(Problem);
      boundary_omega_dwdx_update(Problem);
      poisson_inner_psi_iterator(Problem);
      inner_u_v_compute(Problem);
    }


}
