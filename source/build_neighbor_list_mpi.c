#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "c_extras.h"
#include "mpi_stubs.c"
#include "master_node.c"
#include "work_nodes.c"

/* 
 * Returns the neighbor list and limits for the list with the
 * radius radius
 *
 * arguments:
 *      list:    The location the neighbor lists is written out to
 *      nvlst:   The to which the number of particles in each list is written out to
 *      radius:  The radius of the list (must be one of the radiuses supplied in the
 *                  initialization of the neighbor list builder. Otherwise, behavior will
 *                  be unpredictable)
 *      max:     The max number of particles in each neighbor list for this list
 */
void get_lists_(int* list, int* nvlst, double radius, int max){
    int i,j;
    int list_num = -1;
    for (i = 0; i < NUM_LISTS; ++i){
        if (radius == THRESHOLDS[i]) {
            list_num = i;
            break;
        }
    }
    if (list_num == -1) return;
    for (i = 0; i < n; ++i) {
        nvlst[i] = nvlsts_stored[list_num * n + i];
    }
    for (i = 0; i < n; ++i) {
        for (j = 0; j < nvlst[i]; ++j){
            list[i * max + j] = neighborlist_stored[i * num_neighbors + j];
        }
    }
}

/*
 * Compare function used for qsort
 */
int doublecompare (const void * a, const void * b)
{
    if ( *(double*)a <  *(double*)b ) return -1;
    if ( *(double*)a == *(double*)b ) return 0;
    if ( *(double*)a >  *(double*)b ) return 1;
}
/*
 * supplies the neighbor list builder with required 
 * parameters. Must be called after MPI_Setup, but before initial_nList_build 
 * 
 * arguments:
 *      num_particles: number of particles in the system
 *      comm_sz: Number of nodes. The value that was returned by MPI_Setup
 *      threshold: the distance under which two particles are considered neighbors
 *      move_threshold: the buffer distance for a particle to be considered as having moved
 *      x_width, y_width, z_width: The dimensions of the system
 * 
 * NOTE: threshold, x_width, y_width, z_width should all be in the same units
 */
void setup_nlist_builder_(int num_particles, int comm_sz, int num_lists, double* thresholds, double move_threshold, double x_width, double y_width, double z_width, int maxvlst){
    
    int x, y, z, i;

    num_work_nodes = comm_sz - 1;
    NUM_LISTS = num_lists;
    MOVE_THRESHOLD_2 = move_threshold;
    
    for (i = 0; i < num_lists; ++i){
        THRESHOLDS[i] = thresholds[i];
    }
    
    qsort (THRESHOLDS, num_lists, sizeof(double), doublecompare);

    XDiv = (int)(x_width/THRESHOLDS[num_lists - 1]);
    YDiv = (int)(y_width/THRESHOLDS[num_lists - 1]);
    ZDiv = (int)(z_width/THRESHOLDS[num_lists - 1]);
    DIV_WIDTH_X = x_width/XDiv;
    DIV_WIDTH_Y = y_width/YDiv;
    DIV_WIDTH_Z = z_width/ZDiv;
    
    WIDTH_X = x_width;
    WIDTH_Y = y_width;
    WIDTH_Z = z_width;
    WIDTH_X2 = x_width/2;
    WIDTH_Y2 = y_width/2;
    WIDTH_Z2 = z_width/2;
    
    
    n = num_particles;
    num_neighbors = maxvlst;
    block_size = max(20, 2*n/(XDiv * YDiv * ZDiv));
    
    send_extra_information();
    
    determine_boundries(comm_sz, 0);
    
    int xdim, ydim, zdim, blockdim;
    xdim = sizeof(Particle***) * XDiv;
    ydim = sizeof(Particle**) * YDiv;
    zdim = sizeof(Particle*) * ZDiv;
    blockdim = sizeof(Particle) * block_size;
    
    particles = malloc(xdim);
    for (x = 0; x < XDiv; ++x){
        particles[x] = malloc(ydim);
        for (y = 0; y < YDiv; ++y){
            particles[x][y] = malloc(zdim);
            for (z = 0; z < ZDiv; ++z){
                particles[x][y][z] = malloc(blockdim);
            }
        }
    }    
    offsetx = WIDTH_X/2;
    offsety = WIDTH_Y/2;
    offsetz = WIDTH_Z/2;
    
    nvlsts_stored =  malloc(sizeof(int) * NUM_LISTS * n);
    neighborlist_stored = malloc(sizeof(int) * n * num_neighbors);
}

/*
 *  does a complete neighbor list build from scratch for all particles. 
 *  "setup_nList_builder" must be called before this, and this must be called at least once before
 *  any "subsequent_nList_build".
 *
 *  variables and parameters:
 *      x_coords, y_coords, z_coords: arrays of doubles, of length num_particles, as supplied in
 *          "setup_nList_builder", representing the x, y, and z coordinates of the particles with the corresponding index
 *      neighbors: the location of the int arrays where the neighbor lists will be written out to. 
 *      nvlist: The number of particles in each list
 */
void initial_nlist_build_(double* x_coords, double* y_coords, double* z_coords){
    int particle_counts[XDiv][YDiv][ZDiv];
    
    int x, y, z;
    for (x = 0; x < XDiv; ++x) {
        for (y = 0; y < YDiv; ++y) {
            for (z = 0; z < ZDiv; ++z) {
                particle_counts[x][y][z] = 0;
            }
        }
    }    

    organize_data_initial(x_coords, y_coords, z_coords, n, particle_counts);
    send_instruction(BUILD_NEIGHBOR_LIST);
    distribute_data();
    receive_neighborlists(neighborlist_stored, nvlsts_stored, n);
}




/*
 * does an update of the supplied neighbor list. The only particles 
 * considered will be ones that have moved more than move_threshold as specified in "setup_nList_builder" 
 * Prior to calling subsequent_nList_build, "setup_nList_build" must have been called to specify
 * parameters for the builder, and "initial_nList_build" must have been called to create the initial
 * full build
 *
 *  arguments:
 *      x_coords, y_coords, z_coords: arrays of doubles, of length num_particles, as supplied in
 *          "setup_nList_builder", representing the x, y, and z coordinates of the particles
 *          with the corresponding index
 *      neighbors: the location of the int arrays where the neighbor lists will be written out to. 
 *		nvlist: The number of particles in each list	
 */
void subsequent_nlist_build_(double* x_coords, double* y_coords, double* z_coords){
    pi
    int particle_counts[XDiv][YDiv][ZDiv];
    int x, y, z;
    for (x = 0; x < XDiv; ++x) {
        for (y = 0; y < YDiv; ++y) {
            for (z = 0; z < ZDiv; ++z) {
                particle_counts[x][y][z] = 0;
            }
        }
    }

    organize_data_new(x_coords, y_coords, z_coords, n, particle_counts);
    send_instruction(BUILD_NEIGHBOR_LIST);
    distribute_data();
    receive_neighborlists(neighborlist_stored, nvlsts_stored, n);
}


/*
 * signals the work nodes to shut down and clean up
 *  and frees memory that was being used by the neighborlist builder
 *
 * arguments:
 *      neighbors: the location of the lists built by the neighbor list builder
 */
void cleanup_nlist_builder_(){
    send_instruction(EXIT);
    int x, y, z;
    for (x = 0; x < XDiv; ++x){
        for (y = 0; y < YDiv; ++y){
            for (z = 0; z < ZDiv; ++z){
                free(particles[x][y][z]);
            }
        }
    }
}


/*
 * "start_work_nodes" sets up and starts processes > 0, which will sit and wait for instructions
 * from proc 0
 *
 * variables and parameters:
 *      my_rank: process rank
 *      comm_sz: number of MPI processes
 *      n: number of particles in the system
 */ 
void start_work_nodes_(int my_rank, int comm_sz, int n){
    receive_extra_information();
    
    Neighbor** neighbors;
    Particle* particles[XDiv][YDiv][ZDiv];
    Particle** moved_list[XDiv][YDiv][ZDiv];
    int moved_list_counts[XDiv][YDiv][ZDiv];
    int instruction, j;
    MPI_Datatype* MPI_Particle;

    int first = 1; //Keeps track of whether this is the first iteration of the neighbor list builder
#pragma omp parallel shared (neighbors, particles, moved_list, moved_list_counts, instruction, node_boundries, my_rank)
    {
        int i = 0;
        num_work_nodes = comm_sz - 1;
        initialize_data_structures_work(n, my_rank, comm_sz, &neighbors, particles, moved_list, moved_list_counts);
        if(node_boundries[my_rank-1][0] != -1){
            get_instructions(my_rank, &instruction); 
            while(instruction != EXIT){
                if (instruction == BUILD_NEIGHBOR_LIST) {
                    if (first) {
                        initial_build_work(particles, my_rank, neighbors);
                    }
                    else {
                        subsequent_build_work(particles, my_rank, moved_list_counts, moved_list, neighbors);
                    }
#pragma omp single
                    {
                        first = 0;
                    }
                }
                else if (instruction == PME){
                    grid_mpole_work(particles);
                }
                get_instructions(my_rank, &instruction);
            }
            clean_up_work(particles, moved_list, neighbors);
        }
    }
}








