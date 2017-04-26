#include "cfd.h"

void scalar_rhs(struct _problem* Problem, int i, int j, double du){
    double omega_loc = (*Problem).omega;
    double h_loc = (*Problem).h;
    double domdy,dom2dy2,domdx,dom2dx2 ;

    domdy = (omega_loc[i][j+1]-omega_loc[i][j-1])/(2*h_loc);
    dom2dy2 = (omega_loc[i][j+1]-2*omega_loc[i][j]+omega_loc[i][j-1])/h_loc^2;

    domdx = (omega_loc[i+1][j]-omega_loc[i-1][j])/(2*h_loc);
    dom2dx2 = (omega_loc[i+1][j]-2*omega_loc[i][j]+omega_loc[i-1][j])/h_loc^2;

    du = domdy + dom2dy2 + domdx + dom2dx2;

}

void first_iteration_omega(struct _problem* Problem, double du_old){
    for(int i=0; i<(*Problem).Nx; i++){
        for(int j = 0; j < (*Problem).Ny ; j++){
            du_old=scalar_rhs((*Problem),i,j,du);
            (*Problem).omega[i][j]=(*Problem).omega[i][j]+du_old*(*Problem).dt;
        }
    }
}

void integration_omega(struct _problem* Problem){
    double du,du_old;
    first_iteration_omega(*Problem,du_old);

    for(int i=0; i<(*Problem).Nx; i++){
        for(int j = 0; j < (*Problem).Ny ; j++){
            du=scalar_rhs((*Problem),i,j,du);
            (*Problem).omega[i][j]=(*Problem).omega[i][j]+dt/2*(3*du-du_old);
        }
    }

}
