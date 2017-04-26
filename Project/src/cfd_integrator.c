#include "cfd.h"

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
