#include "cfd.h"
#include <unistd.h>
#include <stdio.h>
#include <stddef.h>

double functionQ(struct _problem* Problem){
  if((*Problem).flag_os == 1)  return (*Problem).Q0*cos(2*M_PI*(*Problem).f*(*Problem).t);
  else                         return (*Problem).Q0;

}

double scalar_rhs_conv(struct _problem* Problem, int i, int j){
    //double**omega_loc = (*Problem).omega;
    double domdx,domdy;

    domdx = (*Problem).u[i][j]*( (*Problem).omega[i+1][j]   -  (*Problem).omega[i-1][j]   )/(2.0*(*Problem).h);
    domdy = (*Problem).v[i][j]*( (*Problem).omega[i]  [j+1] -  (*Problem).omega[i]  [j-1] )/(2.0*(*Problem).h);

    //if( j == (*Problem).NHs + 5 && i == (*Problem).NLs ) fprintf(stderr, "OMEGA %f %f \n", domdy, domdx );
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
   (*Problem).t   = (*Problem).t + (*Problem).dt;
   (*Problem).tau = (*Problem).tau + (*Problem).dtau;
    for(int i=1; i<(*Problem).Nx-1; i++){
        for(int j = (*Problem).imap[i]+1; j < (*Problem).Ny-1; j++){
            dom_old_conv = scalar_rhs_conv(Problem,i,j);
            dom_old_diff = scalar_rhs_diff(Problem,i,j);
            //if( j == (*Problem).NHs + 3 && i == (*Problem).NLs - 2 ) fprintf(stderr, "%f %f \n", dom_old_conv,(*Problem).omega[(*Problem).NLs - 2][(*Problem).NHs + 3] );
            (*Problem).omega[i][j]= (*Problem).omega[i][j] - (*Problem).dt*dom_old_conv + (*Problem).dt*dom_old_diff;
            (*Problem).f_old[i][j]= dom_old_conv;
        }
    }
    boundary_psi_update(Problem,functionQ);
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

    for(int k=0; k < (*Problem).Ntime+1; k++ ){
      (*Problem).t = (k+1)*(*Problem).dt;
      (*Problem).tau = (k+1)*(*Problem).dtau;
      // ## Calcule Omega+1
      for(int i=1; i<(*Problem).Nx-1; i++){
          for(int j = (*Problem).imap[i]+1; j < (*Problem).Ny-1; j++){
              dom_conv=scalar_rhs_conv(Problem,i,j);
              dom_diff=scalar_rhs_diff(Problem,i,j);
              //if( j == (*Problem).NHs + 3 && i == (*Problem).NLs - 2 ) fprintf(stderr, "%f %f \n", dom_conv,(*Problem).omega[(*Problem).NLs - 2][(*Problem).NHs + 3] );
              (*Problem).omega[i][j] = (*Problem).omega[i][j] - (*Problem).dt*0.5*(3.0*dom_conv - (*Problem).f_old[i][j] ) + (*Problem).dt*dom_diff ;
              (*Problem).f_old[i][j] = dom_conv;

              check = diagnose_check(Problem,i,j,k);
              //if(check != 0 ){ fprintf(stderr, "[DEADSTOP] Diagnose Check: %d\n",check); deadstop_exit(Problem);}
          }
      }
      boundary_psi_update(Problem,functionQ);
      poisson_inner_psi_iterator(Problem);
      boundary_omega_update(Problem);
      boundary_omega_dwdx_update(Problem);
      inner_u_v_compute(Problem);

      snapshot(Problem);

    }
}

void snapshot(struct _problem* Problem){
  if( (*Problem).t_snapshot != 0 &&((int)(100*(*Problem).tau/(*Problem).dtau))%( (int) ((*Problem).t_snapshot*100/(*Problem).dtau)) == 0){
    print_problem_data(Problem);
    if((*Problem).flag_pres == 1){
      u_v_stag(Problem);
      boundary_pression_ghost_update(Problem);
      boundary_pression_ghost_in_out(Problem);
      poisson_inner_pres_iterator(Problem);
      print_problem_pressure(Problem);
    }
  }
}
