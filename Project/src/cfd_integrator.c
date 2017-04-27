#include "cfd.h"

void time_integration(struct _problem* Problem){
	(*Problem).t = (*Problem).t + dt;
	printf("omega\n");

}

void first_time_integration(struct _problem* Problem){
	(*Problem).t = (*Problem).t + (*Problem).dt;

	printf("poisson\n");
	poisson_inner_psi_iterator(*Problem);

	first_iteration_omega(*Problem);





	integration_omega((*Problem));
	inner_u_v_compute((*Problem));
}

void firstIteration(problem *myProb)
{


    //firstIterationUstar(myProb);
    //firstIterationVstar(myProb);


    poisson(myProb);


    updateU(myProb);
    updateV(myProb);
    updateP(myProb);
    computew(myProb);
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
    double dom,dom_old;
    first_iteration_omega((*Problem));

    for(int i=0; i<(*Problem).Nx; i++){
        for(int j = 0; j < (*Problem).Ny ; j++){
            dom_conv=scalar_rhs_conv((*Problem),i,j);
            dom_diff=scalar_rhs_diff((*Problem),i,j);
            (*Problem).omega[i][j]=(*Problem).omega[i][j]+(*Problem).dt*(1/2*(3*dom_conv-dom_old)+dom_diff);
        }
    }

}
