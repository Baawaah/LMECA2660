#include "cfd.h"
#include <unistd.h>
#include <stdio.h>
#include <stddef.h>
/*
void time_integration(struct _problem* Problem){
	(*Problem).t = (*Problem).t + dt;
	printf("omega\n");
}*/


double functionQ(struct _problem* Problem){
  return  (*Problem).Q0; //(*Problem).Q0*cos(2*M_PI*t);
}
/*
void first_time_integration(struct _problem* Problem){
	(*Problem).t = (*Problem).t + (*Problem).dt;
	//printf("poisson\n");
  boundary_psi_update(Problem,functionQ);
  boundary_omega_update(Problem);
  first_iteration_omega(Problem);
  poisson_inner_psi_iterator(Problem);
	inner_u_v_compute(Problem);
	//integration_omega((*Problem));
}*/
/*
double scalar_rhs_conv(struct _problem* Problem, int i, int j){
    double**omega_loc = (*Problem).omega;
    double h_loc = (*Problem).h;
    double domdx,domdy;
    domdx = (*Problem).u[i][j]*(omega_loc[i+1][j]-omega_loc[i-1][j])/(2*h_loc);
    domdy = (*Problem).v[i][j]*(omega_loc[i][j+1]-omega_loc[i][j-1])/(2*h_loc);
    return -(domdy + domdx);
}*/
double scalar_rhs_conv(struct _problem* Problem, int i, int j){
    double**omega_loc = (*Problem).omega;
    double h_loc = (*Problem).h;
    double domdx,domdy;

    domdx = (*Problem).u[i][j]*( omega_loc[i+1][j]   -  omega_loc[i-1][j]   )/(2.0*(*Problem).h);
    domdy = (*Problem).v[i][j]*( omega_loc[i]  [j+1] -  omega_loc[i]  [j-1] )/(2.0*(*Problem).h);

    return (domdy + domdx);

}

double scalar_rhs_diff(struct _problem* Problem, int i, int j){
    double**omega_loc = (*Problem).omega;
    double h_loc = (*Problem).h;
    double dom2dx2,dom2dy2 ;

    dom2dx2 = (*Problem).nu*( omega_loc[i+1][j]  - 2.0*omega_loc[i][j] + omega_loc[i-1][j]   )/pow((*Problem).h,2.0);
    dom2dy2 = (*Problem).nu*( omega_loc[i]  [j+1]- 2.0*omega_loc[i][j] + omega_loc[i]  [j-1] )/pow((*Problem).h,2.0);

    return (dom2dy2 + dom2dx2);

}
/*
double scalar_rhs_conv_old(struct _problem* Problem, int i, int j){
    //double** omega_old = (*Problem).w_old;
    double h_loc = (*Problem).h;
    double domdx_old,domdy_old;
    domdx_old = (*Problem).u[i][j]*(omega_old[i+1][j]-omega_old[i-1][j])/(2.0*h_loc);
    domdy_old = (*Problem).v[i][j]*(omega_old[i][j+1]-omega_old[i][j-1])/(2.0*h_loc);
    return (domdy_old + domdx_old);
}*/

// double scalar_rhs_diff_old(struct _problem* Problem, int i, int j){
//     double omega_old = (*Problem).w_old;
//     double h_loc = (*Problem).h;
//     double dom2dx2_old,dom2dy2_old ;

//     dom2dx2_old = (*Problem).nu*(omega_old[i+1][j]-2*omega_old[i][j]+omega_old[i-1][j])/h_loc^2;
//     dom2dy2_old = (*Problem).nu*(omega_old[i][j+1]-2*omega_old[i][j]+omega_old[i][j-1])/h_loc^2;

//     return (dom2dy2_old + dom2dx2_old);

// }

void first_iteration_omega(struct _problem* Problem){
	double dom_old_conv,dom_old_diff;
    for(int i=1; i<(*Problem).Nx-1; i++){
        for(int j = 1; j < (*Problem).imax_map[i]-1; j++){
            dom_old_conv=scalar_rhs_conv(Problem,i,j);
            dom_old_diff=scalar_rhs_diff(Problem,i,j);
            //printf(stderr,"GNA: %f /n",dom_old_conv);
            (*Problem).omega[i][j]=(*Problem).omega[i][j] -(*Problem).dt*dom_old_conv
                                                          +(*Problem).dt*dom_old_diff;
            (*Problem).f_old[i][j]= dom_old_conv;
        }
    }
    boundary_omega_update(Problem);
    boundary_psi_update(Problem,functionQ);
    poisson_inner_psi_iterator(Problem);
    inner_u_v_compute(Problem);
}

void integration_omega(struct _problem* Problem){
    int check = 0;
    double dom_conv,dom_diff;
    //double** w_old = (*Problem).w_old;
    first_iteration_omega(Problem);
    //print_problem_data(Problem);
    for(int k=1; k < (*Problem).Ntime; k++ ){
      (*Problem).t = (*Problem).t + (*Problem).dt;
      // ## Calcule Omega+1
      for(int i=1; i<(*Problem).Nx-1; i++){
          for(int j = 1; j < (*Problem).imax_map[i]-1; j++){
              dom_conv=scalar_rhs_conv(Problem,i,j);
              dom_diff=scalar_rhs_diff(Problem,i,j);
              (*Problem).omega[i][j] = (*Problem).omega[i][j] - (*Problem).dt*0.5*(3.0*dom_conv - (*Problem).f_old[i][j] )
                                                              + (*Problem).dt*dom_diff ;
              (*Problem).f_old[i][j] = dom_conv;

              check = diagnose_check(Problem,i,j);
              //if(check != 0 ){ fprintf(stderr, "[DEADSTOP] Diagnose Check: %d\n",check); deadstop_exit(Problem);}
          }
      }

      poisson_inner_psi_iterator(Problem);

      boundary_omega_update(Problem);

      boundary_psi_update(Problem,functionQ);

      inner_u_v_compute(Problem);
    }


}