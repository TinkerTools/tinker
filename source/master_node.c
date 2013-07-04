#include "stdlib.h"

/*
 * Determines the block boundries for each node, writen out to global
 * node_boundries
 * 
 * comm_sz: Takes in the number of work nodes avaiable, returns
 *                the number of work nodes actually in use
 * my_rank: node rank
 */
void determine_boundries(int comm_sz, int my_rank){
    int i;
    num_work_nodes = comm_sz - 1;

    node_boundries = malloc(sizeof(int*) * num_work_nodes);
    for (i = 0; i < num_work_nodes; ++i) {
        node_boundries[i] = malloc(sizeof(int) * 6);
    }
    
    int full_nodes = num_work_nodes;
    
    if (num_work_nodes > (XDiv * YDiv * ZDiv)){
        num_work_nodes = (XDiv * YDiv * ZDiv);
    }
    
    int base = (int) log2(num_work_nodes);
    if (base != log2(num_work_nodes)){
        num_work_nodes = pow(2, base);
    }
    
    int XSplit = (log2(num_work_nodes) + 2) / 3;
    int YSplit = (log2(num_work_nodes) + 1) / 3;
    int ZSplit = (log2(num_work_nodes) + 0) / 3;  
    
    int X2 = (int) (pow(2, XSplit));
    int Y2 = (int) (pow(2, YSplit));
    int Z2 = (int) (pow(2, ZSplit));
    
    int XLength = XDiv / X2; 
    int YLength = YDiv / Y2;
    int ZLength = ZDiv / Z2;    
    
    
    for (i = 0; i < num_work_nodes; ++i) {
        node_boundries[i][0] = (i % X2)  * XLength; 
        
        if ((i % X2) == X2 - 1) node_boundries[i][1] = XDiv - 1;
        else node_boundries[i][1] = node_boundries[i][0] + XLength - 1; 
        
        node_boundries[i][2] = ((i/X2) % Y2) * YLength; 
        
        if ((i/X2) % Y2 == Y2 - 1) node_boundries[i][3] = YDiv - 1;
        else node_boundries[i][3] = node_boundries[i][2] + YLength - 1;  
        
        node_boundries[i][4] = ((i/(X2 * Y2)) % Z2) * ZLength;  
        
        if ((i/(X2 * Y2)) % Z2 == Z2 - 1) node_boundries[i][5] = ZDiv - 1;
        else node_boundries[i][5] = node_boundries[i][4] + ZLength - 1; 
        
    }

    for (i = num_work_nodes; i < full_nodes; ++i) {
        node_boundries[i][0] = -1;
    }    
}

/*
 * Spatially decomposes the particle location information into discrete blocks, to prepare for sending
 * to the work nodes
 *
 * arguments:
 *      x_coords_new: The x coordinates of the particles
 *      y_coords_new: The y coordinates of the particles
 *      z_coords_new: The z coordinates of the particles
 *      n:            The number of particles
 *      indexes(OUT): The number of particles in each block
 */
void organize_data_initial(double* x_coords_new, double* y_coords_new, double* z_coords_new, int n, int indexes[XDiv][YDiv][ZDiv]){
    
    int x, y, z, i;
    Particle* p;
    
    // Indexes keeps track of how many particles are currently in each block,
    // That is, the next location to put a new particle in in each block
    for (x = 0; x < XDiv; ++x) {
        for (y = 0; y < YDiv; ++y) {
            for (z = 0; z < ZDiv; ++z){
                indexes[x][y][z] = 0;
            }
        }
    }
    
    //
    for (i = 0; i < n; ++i){                       
        image_(&x_coords_new[i], &y_coords_new[i], &z_coords_new[i]);
        
        //Calculate which block each particle belongs in based on location
        //Since the particles can be anywhere in space, but the blocks start
        //with 0,0,0, the offsets move the particles, such that they all have
        //postive coordinates and the box starts at 0,0,0
        x = (x_coords_new[i] + offsetx)/WIDTH_X * XDiv;
        y = (y_coords_new[i] + offsety)/WIDTH_Y * YDiv;
        z = (z_coords_new[i] + offsetz)/WIDTH_Z * ZDiv; 

        p =  &particles[x][y][z][indexes[x][y][z]];
        
        p->pos[0] = x_coords_new[i];
        p->pos[1] = y_coords_new[i];
        p->pos[2] = z_coords_new[i];
        
        p->oldPos[0] = x_coords_new[i];
        p->oldPos[1] = y_coords_new[i];
        p->oldPos[2] = z_coords_new[i];

        indexes[x][y][z]++;
        
        p->index = invIndex(i);
        p->moved = 0;
    }

    // Mark that all the empty spots are not valid particles
    for (x = 0; x < XDiv; ++x) {
        for (y = 0; y < YDiv; ++y) {
            for (z = 0; z < ZDiv; ++z){
                for (i = indexes[x][y][z]; i < block_size; ++i) {
                    particles[x][y][z][i].index = -1;
                }
            }
        }
    }
}

/*
 * Organize the particle location information in every iteration past the initial
 *
 * arguments:
 *      x_coords_new: The x coordinates of the particles
 *      y_coords_new: The y coordinates of the particles
 *      z_coords_new: The z coordinates of the particles
 *      n:            The number of particles
 *      indexes(OUT): The number of particles in each block
 */
void organize_data_new(double* x_coords_new, double* y_coords_new, double* z_coords_new, int n, int indexes[XDiv][YDiv][ZDiv]){
    //Iterate through blocks, update each location. Move between blocks if necessary
    
    int i, j, x, y, z;
    int xnew, ynew, znew;
    int index;
    Particle* p;
    
    // Mark all particles as not having moved to reset the previous iteration's information
    for (x = 0; x < XDiv; ++x) {
        for (y = 0; y < YDiv; ++y) {
            for (z = 0; z < ZDiv; ++z){
                for (i = 0; i < block_size; ++i) {
                    particles[x][y][z][i].moved = 0;
                }
            }
        }
    }
    
    // Instead of placing particles in blocks, this goes through each particle already in a block
    // updates its location information, and moves it between blocks if necessary
    for (x = 0; x < XDiv; ++x) {
        for (y = 0; y < YDiv; ++y) {
            for (z = 0; z < ZDiv; ++z) {
                for (i = 0; i < block_size; ++i) {
                    
                    p = &particles[x][y][z][i];
                    if (p->index != -1 && p->moved != 1){

                        index = pIndex(p->index);
                        image_(&x_coords_new[index], &y_coords_new[index], &z_coords_new[index]);
                        xnew = (x_coords_new[index] + offsetx)/WIDTH_X * XDiv;
                        ynew = (y_coords_new[index] + offsety)/WIDTH_Y * YDiv;
                        znew = (z_coords_new[index] + offsetz)/WIDTH_Z * ZDiv;
                        

                        if (x != xnew || y != ynew || z != znew){
                            //Find a new empty location for the particle
                            for (j = 0; j < block_size; ++j) {
                                if (particles[xnew][ynew][znew][j].index == -1) break; 
                            }
                            particles[xnew][ynew][znew][j] = particles[x][y][z][i];
                            particles[xnew][ynew][znew][j].moved = 1; 
                            particles[x][y][z][i].index = -1;
                            p = &particles[xnew][ynew][znew][j];

                        }
                        p->pos[0] = x_coords_new[index];
                        p->pos[1] = y_coords_new[index];
                        p->pos[2] = z_coords_new[index];
                        double xDif = p->pos[0] - p->oldPos[0]; //The difference between the current
                        double yDif = p->pos[1] - p->oldPos[1]; //position and the position the last 
                        double zDif = p->pos[2] - p->oldPos[2]; //time that particle's list was rebuilt
                        
                        if ((xDif*xDif) + (yDif* yDif) + (zDif * zDif) > MOVE_THRESHOLD_2) {
                            p->moved = 1;
                        }
                        if (p->moved == 1) {
                            p->oldPos[0] = x_coords_new[index];  
                            p->oldPos[1] = y_coords_new[index];
                            p->oldPos[2] = z_coords_new[index];
                        }                        
                    }    
                }
            }
        }
    }
}

/*
 * Distributes updated particle data to all work nodes
 */
void distribute_data(){
    int i, x, y, z;
    int XNew, YNew, ZNew;
    for (i = 0; i < num_work_nodes; ++i) {
        for (x = node_boundries[i][0]; x <= node_boundries[i][1]+1; ++x) {
            for (y = node_boundries[i][2]-1; y <= node_boundries[i][3]+1; ++y) {
                for (z = node_boundries[i][4]-1; z <= node_boundries[i][5]+1; ++z) {
                    XNew = (x+XDiv)%XDiv;
                    YNew = (y+YDiv)%YDiv;
                    ZNew = (z+ZDiv)%ZDiv;
                    int indexes = XNew+100*YNew+10000*ZNew; //A unique tag for each block
                    MPI_Send(particles[XNew][YNew][ZNew], block_size, MPI_Particle, i+1, indexes, MPI_COMM_WORLD);                        
                    
                }   
            }
        }
    }

}

/*
 * Receives the fully built neighbor lists from work nodes and decompresses the message
 *
 * arguments:
 *      neighbors:  The location that the completed neighbor lists will be stored
 *      nvlst:      The number of particles in each list
 *      n:          The total number of particles in the system
 */
void receive_neighborlists(int* neighbors, int* nvlst, int n){
    int i, j, k, x, y, z, tag;
    int index, array_index;
    
    int size = num_neighbors * (n / num_work_nodes);
    int* message = malloc(sizeof(int) * size);
    int* nvlst_temp = malloc(sizeof(int) * n * NUM_LISTS);
    for (i = 0; i < num_work_nodes; ++i){
        MPI_Recv(message, size, MPI_INT, i+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(nvlst_temp, NUM_LISTS * n, MPI_INT, i+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        k = 0;
        
        while (k < size) {
            // End of Message
            if (message[k] == -2) {
                break;   
            }
            // End of single neighborlist
            else if (message[k] == -1) {
                ++k;
                if (k < size) index = message[k];
                ++k;
            }
            
            // Copy nvlst for the particles each work node was responsible for
            for (j = 0; j < NUM_LISTS; ++j){
                nvlst[j*n + pIndex(index)] = nvlst_temp[j*n + pIndex(index)];
            }
            
            //Copy the neighbor list
            j = 0;
            while(message[k] != -1 && message[k] != -2 && j < num_neighbors && k < size){
                if (j >= num_neighbors) printf("Exceeded size of neighborlist. Increse maxvlst");
                array_index = pIndex(index) * num_neighbors + j ;
                neighbors[array_index] = message[k];
                ++k; ++j;
            }
        }
    }
    free(message);
    free(nvlst_temp);
}


/*
 * Send an instruction to the work nodes. Instruction codes are defined in c_extras.h
 */
void send_instruction(int instruction){
    int i;
    for (i = 0; i < num_work_nodes; ++i){
        MPI_Send(&instruction, 1, MPI_INT, i+1, 0, MPI_COMM_WORLD);
    }

}

/*
 * Send various parameters needed for the work nodes
 */
void send_extra_information(){
    int i;
    int int_info[6] = {XDiv, YDiv, ZDiv, n, num_neighbors, NUM_LISTS};
    double double_info[3] = {DIV_WIDTH_X, DIV_WIDTH_Y, DIV_WIDTH_Z};
    
    for (i = 1; i <= num_work_nodes; ++i){
        MPI_Send(int_info, 6, MPI_INT, i, 0, MPI_COMM_WORLD);
        MPI_Send(double_info, 3, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        MPI_Send(THRESHOLDS, 10, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    }
}