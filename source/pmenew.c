#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "c_extras.h"
int num_work_nodes;
void send_instruction(int instruction);

//The following calculate the index of an item from the fortran array indexes
//a,b... are the indexes and limA, limB... are the dimensions of the array
#define index2(a,b,limA,limB) ((b)*limA + (a))
#define index3(a,b,c,limA,limB,limC) ((c)*limA*limB + (b)*limA + (a))
#define index4(a,b,c,d,limA,limB,limC,limD) ((d)*limA*limB*limC + (c)*limA*limB + (b)*limA + (a))

int bsorder, nfft1, nfft2, nfft3, npole;

/*
 *  The c equivalent of the isign function
 */
int intsign(int a, int b){
    if (a * b >= 0) return a;
    else return -a;
}

/*
 * Performs the calculations for a single particle for grid_mpole
 * Equivalent to the calculations in the grid_mpole method in pmestuff.f
 */
void mpole_math(double* fmp, int* igrid, double* qgrid, double* thetai1, double* thetai2, double* thetai3, int bsorder, int nfft1, int nfft2, int nfft3, int npole, int n, int m){
    int i, j, k;
    int i0, j0, k0;
    int it1, it2, it3;
    int igrd0, jgrd0, kgrd0;
    double v0, u0, t0;
    double v1, u1, t1;
    double v2, u2, t2;
    double term0, term1, term2;
    
    igrd0 = igrid[index2(0, m, 3, maxatm)];
    jgrd0 = igrid[index2(1, m, 3, maxatm)];
    kgrd0 = igrid[index2(2, m, 3, maxatm)];
    k0 = kgrd0;
    
    for (it3 = 0; it3 < bsorder; ++it3){
        k0 = k0 + 1;
        k = k0 + 1 + (nfft3-intsign(nfft3,k0))/2;
        v0 = thetai3[index3(0,it3,m, 4,bsorder,n)];
        v1 = thetai3[index3(1,it3,m, 4,bsorder,n)];
        v2 = thetai3[index3(2,it3,m, 4,bsorder,n)];
        j0 = jgrd0;
        
        for (it2 = 0; it2 < bsorder; ++it2){
            j0 = j0 + 1;
            j = j0 + 1 + (nfft2-intsign(nfft2,j0))/2;
            u0 = thetai2[index3(0,it2,m, 4,bsorder,n)];
            u1 = thetai2[index3(1,it2,m, 4,bsorder,n)];
            u2 = thetai2[index3(2,it2,m, 4,bsorder,n)];
            term0 = fmp[index2(0,m, 10, maxatm)]*u0*v0 + fmp[index2(2,m, 10, maxatm)]*u1*v0 + fmp[index2(3,m, 10, maxatm)]*u0*v1 + fmp[index2(5,m, 10, maxatm)]*u2*v0 + fmp[index2(6,m, 10, maxatm)]*u0*v2 + fmp[index2(9,m, 10, maxatm)]*u1*v1;
            term1 = fmp[index2(1,m, 10, maxatm)]*u0*v0 + fmp[index2(7,m, 10, maxatm)]*u1*v0 + fmp[index2(8,m, 10, maxatm)]*u0*v1;
            term2 = fmp[index2(4,m, 10, maxatm)] * u0 * v0;
            i0 = igrd0;
            
            for (it1 = 0; it1 < bsorder; ++it1){
                i0 = i0 + 1;
                i = i0 + 1 + (nfft1-intsign(nfft1,i0))/2;
                
                t0 = thetai1[index3(0,it1,m, 4,bsorder,n)];
                t1 = thetai1[index3(1,it1,m, 4,bsorder,n)];
                t2 = thetai1[index3(2,it1,m, 4,bsorder,n)];
                
                qgrid[index4(0,i-1,j-1,k-1,2, nfft3,nfft2,nfft1)] = qgrid[index4(0,i-1,j-1,k-1,2, nfft3,nfft2,nfft1)] + term0*t0 + term1*t1 + term2*t2;
            }
        }
    }
    
}

/*
 * Master node sends information that doesn't change between
 */
void send_single_info(int bsorder, int nfft1, int nfft2, int nfft3, int npole, int n){
    int i;
    int info[5] = {bsorder, nfft1, nfft2, nfft3, npole};//, n};
    for (i = 1; i <= num_work_nodes; ++i){
        MPI_Send(info, 5, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
    
}

/*
 * Work node recieves the information sent by "send_single_info"
 */
void recieve_single_info(){
    int info[5];
    MPI_Recv(info, 5, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    bsorder = info[0];
    nfft1   = info[1];
    nfft2   = info[2];
    nfft3   = info[3];
    npole   = info[4];
}


/*
 * Master node sends information needed to perform a single iteration
 */
void send_step_info(double* fmp, int* igrid, double* thetai1, double* thetai2, double* thetai3, int npole, int n, int bsorder){
    int i;
    for (i = 1; i <= num_work_nodes; ++i){
        MPI_Send(fmp, 10 * npole, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        MPI_Send(igrid, 10 * n, MPI_INT, i, 0, MPI_COMM_WORLD);
        MPI_Send(thetai1, 4 * bsorder * n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        MPI_Send(thetai2, 4 * bsorder * n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        MPI_Send(thetai3, 4 * bsorder * n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    }
}

/*
 * Work node recieves information sent by "send_step_info"
 */
void recieve_step_info(double* fmp, int* igrid, double* thetai1, double* thetai2, double* thetai3){
    MPI_Recv(fmp, 10 * npole, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(igrid, 10 * n, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(thetai1, 4 * bsorder * n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(thetai2, 4 * bsorder * n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(thetai3, 4 * bsorder * n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

/*
 * Work node sends its portion of the completed grid back to the master node
 */
void send_grid(double* qgrid){
    int size = 2*nfft1*nfft2*nfft3;
    MPI_Send(qgrid, size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
}

/*
 * Master node recieves and reduces completed grids
 */
void receive_grids(double* qgrid, int nfft1, int nfft2, int nfft3){

    int i, j;
    int size = 2*nfft1*nfft2*nfft3;
    double* tempgrid = malloc(size * sizeof(double));
    for (i = 1; i <= num_work_nodes; ++i){
        MPI_Recv(tempgrid, size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (j = 0; j < size; ++j){
            qgrid[i] += tempgrid[i];
        }
    }
    
    free(tempgrid);
}


int first_pme = 1;
void grid_mpole_mpi_(double* fmp, int* igrid, double* qgrid, double* thetai1, double* thetai2, double* thetai3, int bsorder, int nfft1, int nfft2, int nfft3, int npole, int n){
    int i;
    int size = 2*nfft1*nfft2*nfft3;
    for (i = 0; i < size; ++i){
        qgrid[i] = 0;
    }

    send_instruction(PME);

    if (first_pme){
        send_single_info(bsorder, nfft1, nfft2, nfft3, npole, n);
        first_pme = 0;
    }
    send_step_info(fmp, igrid, thetai1, thetai2, thetai3, npole, n, bsorder);
    receive_grids(qgrid, nfft1, nfft2, nfft3);
}

void grid_mpole_work(Particle* particles[XDiv][YDiv][ZDiv]){

    double *fmp, *thetai1, *thetai2, *thetai3;
    int *igrid;

    first_pme = 1;
    if (first_pme){
        recieve_single_info();
        first_pme = 0;
    }
    
    fmp = malloc(10 * npole * sizeof(double));
    igrid = malloc(10 * n * sizeof(int));
    thetai1 = malloc(4 * bsorder * n * sizeof(double));
    thetai2 = malloc(4 * bsorder * n * sizeof(double));
    thetai3 = malloc(4 * bsorder * n * sizeof(double));
    
    double* qgrid = malloc(2 * nfft1 * nfft2 * nfft3 * sizeof(double));
    
    recieve_step_info(fmp, igrid, thetai1, thetai2, thetai3);

    int x, y, z, i;
    for (x = node_boundries[my_rank][0]; x <= node_boundries[my_rank][1]; ++x) {
        for (y = node_boundries[my_rank][2]; y <= node_boundries[my_rank][3]; ++y) {
            for (z = node_boundries[my_rank][4]; z <= node_boundries[my_rank][5]; ++z) {
                for (i = 0; i < block_size; ++i){
                    mpole_math(fmp, igrid, qgrid, thetai1, thetai2, thetai3, bsorder, nfft1, nfft2, nfft3, npole, n, particles[x][y][z][i].index);
                }
            }
        }
    }
    send_grid(qgrid);
    free(fmp); free(igrid); free(thetai1); free(thetai2); free(thetai3); free(qgrid);
}

