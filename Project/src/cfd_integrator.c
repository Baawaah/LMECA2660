#include "cfd.h"

double scalar_rhs_conv(struct _problem* Problem, int i, int j){
    double omega_loc = (*Problem).omega;
    double h_loc = (*Problem).h;
    double domdx,domdy;

    domdx = (omega_loc[i+1][j]-omega_loc[i-1][j])/(2*h_loc);
    domdy = (omega_loc[i][j+1]-omega_loc[i][j-1])/(2*h_loc);

    return (domdy + domdx);

}

double scalar_rhs_diff(struct _problem* Problem, int i, int j){
    double omega_loc = (*Problem).omega;
    double h_loc = (*Problem).h;
    double dom2dx2,dom2dy2 ;

    dom2dx2 = (omega_loc[i+1][j]-2*omega_loc[i][j]+omega_loc[i-1][j])/h_loc^2;
    dom2dy2 = (omega_loc[i][j+1]-2*omega_loc[i][j]+omega_loc[i][j-1])/h_loc^2;

    return (dom2dy2 + dom2dx2);

}

void first_iteration_omega(struct _problem* Problem, double du_old){
    for(int i=0; i<(*Problem).Nx; i++){
        for(int j = 0; j < (*Problem).Ny ; j++){
            du_old_conv=scalar_rhs_conv((*Problem),i,j);
            du_old_diff=scalar_rhs_diff((*Problem),i,j);
            (*Problem).omega[i][j]=(*Problem).omega[i][j]+(*Problem).dt*(du_old_conv+du_old_diff); // EE for both at the first time step
        }
    }
}

void integration_omega(struct _problem* Problem){
    double du,du_old;
    first_iteration_omega(*Problem,du_old);

    for(int i=0; i<(*Problem).Nx; i++){
        for(int j = 0; j < (*Problem).Ny ; j++){
            du_conv=scalar_rhs_conv((*Problem),i,j);
            du_diff=scalar_rhs_diff((*Problem),i,j);
            (*Problem).omega[i][j]=(*Problem).omega[i][j]+(*Problem).dt*(1/2*(3*du_conv-du_old)+du_diff);
        }
    }

}
