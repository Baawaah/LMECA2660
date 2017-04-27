#include "cfd.h"

void time_integration(struct _problem* Problem){
	(*Problem).t = (*Problem).t + dt;
	printf("omega\n");

}

void first_time_integration(struct _problem* Problem){
	(*Problem).t = (*Problem).t + (*Problem).dt;
	
	printf("u and v\n");
	inner_u_v_compute((*Problem));

	printf("poisson\n");
	poisson_inner_psi_iterator((*Problem));

	first_iteration_omega((*Problem));
	inner_u_v_old_compute((*Problem));

	//integration_omega((*Problem));
}

double scalar_rhs_conv(struct _problem* Problem, int i, int j){
    double omega_loc = (*Problem).omega;
    double h_loc = (*Problem).h;
    double domdx,domdy;

    domdx = (*Problem).u[i][j]*(omega_loc[i+1][j]-omega_loc[i-1][j])/(2*h_loc);
    domdy = (*Problem).v[i][j]*(omega_loc[i][j+1]-omega_loc[i][j-1])/(2*h_loc);

    return -(domdy + domdx);

}

double scalar_rhs_diff(struct _problem* Problem, int i, int j){
    double omega_loc = (*Problem).omega;
    double h_loc = (*Problem).h;
    double dom2dx2,dom2dy2 ;

    dom2dx2 = (*Problem).nu*(omega_loc[i+1][j]-2*omega_loc[i][j]+omega_loc[i-1][j])/h_loc^2;
    dom2dy2 = (*Problem).nu*(omega_loc[i][j+1]-2*omega_loc[i][j]+omega_loc[i][j-1])/h_loc^2;

    return (dom2dy2 + dom2dx2);

}

double scalar_rhs_conv_old(struct _problem* Problem, int i, int j){
    double omega_old = (*Problem).w_old;
    double h_loc = (*Problem).h;
    double domdx_old,domdy_old;

    domdx_old = (*Problem).u[i][j]*(omega_old[i+1][j]-omega_old[i-1][j])/(2*h_loc);
    domdy_old = (*Problem).v[i][j]*(omega_old[i][j+1]-omega_old[i][j-1])/(2*h_loc);

    return -(domdy_old + domdx_old);

}

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
    for(int i=0; i<(*Problem).Nx; i++){
        for(int j = 0; j < (*Problem).Ny ; j++){
            dom_old_conv=scalar_rhs_conv((*Problem),i,j);
            dom_old_diff=scalar_rhs_diff((*Problem),i,j);
            //(*Problem).omega[i][j]=(*Problem).omega[i][j]+(*Problem).dt*(dom_old_conv+dom_old_diff); // EE for both at the first time step
            (*Problem).w_old[i][j]=(*Problem).omega[i][j]+(*Problem).dt*(dom_old_conv+dom_old_diff); 
        }
    }
}

void integration_omega(struct _problem* Problem){
    double dom_conv,dom_diff;
    double w_old = (*Problem).w_old;
    first_iteration_omega((*Problem));

    for(int i=0; i<(*Problem).Nx; i++){
        for(int j = 0; j < (*Problem).Ny ; j++){
            dom_conv=scalar_rhs_conv_old((*Problem),i,j);
            dom_diff=scalar_rhs_diff((*Problem),i,j);
            (*Problem).omega[i][j]=(*Problem).omega[i][j]+(*Problem).dt*(1/2*(3*dom_conv-w_old[i][j])+dom_diff);
        }
    }

}
