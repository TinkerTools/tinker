#include <mpi.h>

void create_type_MPI_PARTICLE(MPI_Datatype* datatype){
	int array_of_block_lengths[4] = {3, 3, 1, 1};
	MPI_Aint a_addr, b_addr, c_addr, d_addr;
	Particle p;
	MPI_Get_address(&p.pos, &a_addr);
	MPI_Get_address(&p.oldPos, &b_addr);
	MPI_Get_address(&p.index, &c_addr);
	MPI_Get_address(&p.moved, &d_addr);
	MPI_Aint array_of_displacements[4] = { 0, b_addr - a_addr, c_addr - a_addr, d_addr - a_addr};
	MPI_Datatype array_of_types[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT};
	MPI_Type_create_struct(4, array_of_block_lengths, array_of_displacements, array_of_types, datatype);
	MPI_Type_commit(datatype);
    
    MPI_Aint size;
    
}


/* "MPI_Setup" initializes MPI and sets up needed MPI Datatypes
 * for communication
 *
 * variables and parameters:
 * OUT my_rank: process rank
 * OUT comm_sz: number of MPI_Processes
 */
void mpi_setup_(int* my_rank, int* comm_sz){
	MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, comm_sz); 
	
	create_type_MPI_PARTICLE(&MPI_Particle);
}

/*
 * "MPI_Cleanup" frees MPI Datatypes and cleans up MPI 
 */ 
void mpi_cleanup_(){
	MPI_Type_free(&MPI_Particle);
	MPI_Finalize();
}

void mpi_cleanup_work_(){
    MPI_Finalize();
}
