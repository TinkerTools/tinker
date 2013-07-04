
/*
 * Clone of the imagen function in nblist.f, for orthagonal systems
 */
void imagen_orthagonal(double* x, double* y, double* z){
    *x = abs(*x);
    *y = abs(*y);
    *z = abs(*z);

    if (*x > WIDTH_X2) *x = *x - WIDTH_X;
    if (*y > WIDTH_Y2) *y = *y - WIDTH_Y;
    if (*z > WIDTH_Z2) *z = *z - WIDTH_Z;

}


/*
 * Checks the distance between two particles, and will return
 * the smallest threshold number that that distance falls under, or 0
 * otherwise
 */
int within_range(Particle one, Particle two){
    
    // Calculate the pairwise distance between coordinates
    // of two particles and take the image of those values
    double difx = one.pos[0] - two.pos[0];
    double dify = one.pos[1] - two.pos[1];
    double difz = one.pos[2] - two.pos[2];
    imagen_orthagonal(&difx, &dify, &difz);
    
    // Calculates the squared distance
    double dist = (difx*difx) + (dify*dify) + (difz*difz);
    
    // Compare the squared distance against the squared
    // thresholds and return the smallest one
    int i;
    for (i = 0; i < NUM_LISTS; ++i){
        if (dist < THRESHOLDS_2[i]){
            return i+1;
        }
    }
    return 0;
}


/*
 * Looks at all the particles within a neighboring block, x,y,z and checks the distance against
 * given particle p
 *
 * arguments:
 *      particles:      The spatially decomposed particle information
 *      x, y, z:        The indexes of the block in the particles array
 *      neighbor_list:  Particle p's neighborlist
 *      p:              The particle whose neighborlist is being built
 *      index:          The next empty spot in the neighborlist
 */
void check_particles(Particle* particles[XDiv][YDiv][ZDiv], int x, int y, int z, Neighbor* neighbor_list, Particle p, int* index) {
    
    // Wraps around to check for "neighboring" particles on the opposite side of the block
    x = (x+XDiv)%XDiv;
    y = (y+YDiv)%YDiv;
    z = (z+ZDiv)%ZDiv;
        
    int i, dist;
    Particle* list = particles[x][y][z];
    
    // Iterate through all the particles in the block to look for neighbors
    for (i = 0; i < block_size; ++i){
        if (list[i].index != -1) { //If valid particles
            dist = within_range(p, list[i]); 
            if (dist != 0){
                // Adds the neighbor to the next available index in the neighborlist
                neighbor_list[*index].ptr = &list[i];
                neighbor_list[*index].index = list[i].index;
                neighbor_list[*index].list = dist;
                *index = *index + 1;
            }
        }
    }
}

/*
 * Looks at all the particles within block x,y,z containing particle p and checks the distance
 *
 * arguments:
 *      particles:      The spatially decomposed particle information
 *      x, y, z:        The indexes of the block in the particles array
 *      list_index:     The index of particle p in the block, x,y,z
 *      neighbor_list:  Particle p's neighborlist
 *      p:              The particle whose neighborlist is being built
 *      index:          The next empty spot in the neighborlist
 */
void check_particles_center(Particle* particles[XDiv][YDiv][ZDiv], int x, int y, int z, int list_index, Neighbor* neighbor_list, Particle p, int *index) {
    
    int i, dist;
    Particle* list = particles[x][y][z];
    // Iterate through all the particles in the block to look for neighbors    
    for (i = list_index + 1; i < block_size; ++i){
        if (list[i].index != -1) {
            dist = within_range(p, list[i]);
            if (dist != 0){
                // Adds the neighbor to the next available index in the neighborlist
                neighbor_list[*index].ptr = &list[i];
                neighbor_list[*index].index = list[i].index;
                neighbor_list[*index].list = dist;
                *index = *index + 1;
            }
        }
    }
}

/*
 * Looks at all the particles in block x,y,z that have moved significantly in the last iteration and updates
 * particle p's neighborlist
 *
 * arguments:
 *      particles:      The spatially decomposed particle information
 *      x, y, z:        The indexes of the block in the particles array
 *      neighbor_list:  Particle p's neighborlist
 *      p:              The particle whose neighborlist is being built
 *      moved_list:     The list of all the particles that have moved significantly
 *      counts:         The number of particles in each moved list
 */
void check_moved_particles(Particle* particles[XDiv][YDiv][ZDiv], int x, int y, int z, Neighbor* neighbor_list, Particle p, Particle** moved_list[XDiv][YDiv][ZDiv], int counts[XDiv][YDiv][ZDiv]) {
    
        // Wraps around to check for "neighboring" particles on the opposite side of the block
        x = (x+XDiv)%XDiv;
        y = (y+YDiv)%YDiv;
        z = (z+ZDiv)%ZDiv;
    
        int i, dist;
        int index = 0;
        Particle** list = moved_list[x][y][z];
    
        // Iterate through all the particles in the moved list to look for neighbors   
        for (i = 0; i < counts[x][y][z]; ++i){
            if (list[i]->index != -1) {
                dist = within_range(p, *list[i]);
                if (dist != 0){
                    //Iterate until we find an empty position in the neighborlist
                    // NOTE: We can't employ the same strategy of just keeping track of the next empty spot
                    // as we did above because we'll have gaps in the neighborlist
                    while(neighbor_list[index].index != -1){
                        index++;
                    }
                    neighbor_list[index].ptr = list[i];
                    neighbor_list[index].index = list[i]->index;
                    neighbor_list[index].list = dist;
                }
            }
        }
}


/*
 * Checks whether one pointer location is larger than a second one
 * Used to check whether one particle is later in an array than the other
 */
int higher_index(Particle *p, Particle *q){ 
    return p - q > 0 ? 1 : 0;
}

/*
 * Looks at all the particles in block x,y,z that have moved significantly in the last iteration and updates
 * particle p's neighborlist
 *
 * arguments:
 *      particles:      The spatially decomposed particle information
 *      x, y, z:        The indexes of the block in the particles array
 *      neighbor_list:  Particle p's neighborlist
 *      p:              The particle whose neighborlist is being built
 *      moved_list:     The list of all the particles that have moved significantly
 *      counts:         The number of particles in each moved list
 */
void check_moved_particles_center(Particle* particles[XDiv][YDiv][ZDiv], int x, int y, int z, Neighbor* neighbor_list, Particle *p, Particle** moved_list[XDiv][YDiv][ZDiv], int counts[XDiv][YDiv][ZDiv]) {
    //Iterate through all moved particles in list
    //Check if the particle is within the radius
    //Add to neighborlist
    
    int i; 
    int index = 0;
    
    Particle** list = moved_list[x][y][z];
    int dist;
    for (i = 0; i < counts[x][y][z]; ++i){ 
        if (list[i]->index != -1) {
            dist = within_range(*p, *list[i]);
            //Since this is for the block that the particle itself is in, it should only add particles that
            //are in array locations that are after the particle, to avoid double counting
            if (higher_index(list[i], p) && dist != 0){

                //Iterate until we find an empty position in the neighborlist
                while(neighbor_list[index].index != -1){
                    index++;
                }
                neighbor_list[index].ptr = list[i];
                neighbor_list[index].index = list[i]->index;
                neighbor_list[index].list = dist;
            }
        }
    }
}

/*
 * Basic radix sort of a neighbor_list, by the list attribute of the neighbor struct
 *
 * arguments:
 *      neighbor_list: The list to be sorted
 */
void sort_list(Neighbor* neighbor_list){

    int indexes[NUM_LISTS]; 
    int list_end[NUM_LISTS];
    Neighbor temp_list[num_neighbors];
    int i;

    // Copy the neighborlist into a temporary array
    for (i = 0; i < num_neighbors; ++i){
        temp_list[i] = neighbor_list[i];
        neighbor_list[i].index = -1;
    }
    // Initialize the temporary count arrays to 0
    for (i = 0; i < NUM_LISTS; ++i) {
        indexes[i] = 0;
        list_end[i] = 0;
    }
    
    // Count the values in each list
    for (i = 0; i < num_neighbors; ++i){
        if (temp_list[i].index != -1){
            list_end[temp_list[i].list - 1] += 1;
        }
    }
 
 
    // Initialize indexes to where each index would start,
    // that is, the cumulative sum of all values before it
    for (i = 1; i < NUM_LISTS; ++i){
        indexes[i] = indexes[i-1] + list_end[i-1];
    }

    // Iterate through the temporary, and using the indexes array to know where
    // to write the next item for each list, write out the sorted list to the neighbor_list
    // array
    for (i = 0; i < num_neighbors; ++i){
        if (temp_list[i].index != -1) {
            neighbor_list[indexes[temp_list[i].list-1]] = temp_list[i];
            indexes[temp_list[i].list-1] += 1;
        }
    }
}


/*
 * Builds the full neighborlist list from scratch for a particle
 * This is used for the initial build for all the particles, then for rebuilds
 * of moved particles
 *
 * arguments:
 *      particles:      The spatially decomposed particle information
 *      neighbor_list:  The particle's neighborlist
 *      x,y,z:          The indecies for the block which the particle belongs to
 *      array_index:    The index of the particle in the list
 */
void build_list_full(Particle* particles[XDiv][YDiv][ZDiv], Neighbor* neighbor_list, int x, int y, int z, int array_index){
    // Iterates through the buckets in in such a way that for every pair of particles which are
    // potential neighbors, one and only one of them will have to check the distance for that pair. This is done by,
    // if we are looking at particle p, only looking at particles further down in the array in the block containing p,
    // and only looking in a "forward" direction (that is, for each pair of blocks that are "across" from eachother around
    // the block containing p, we consider only of those blocks)

    
    int j, k;
    int index = 0;
    Particle p = particles[x][y][z][array_index];
    
    for (j = -1; j <= 1; ++j){ //y offset
        for (k = -1; k <= 1; ++k) { //z offset
            check_particles(particles, x+1, y+j, z+k, neighbor_list, p, &index);
        }
    }
    
    check_particles(particles, x, y+1, z+1, neighbor_list, p, &index);
    check_particles(particles, x, y+1, z  , neighbor_list, p, &index);
    check_particles(particles, x, y+1, z-1, neighbor_list, p, &index);
    check_particles(particles, x, y  , z+1, neighbor_list, p, &index);
    check_particles_center(particles, x, y, z, array_index, neighbor_list, p, &index);
    
    sort_list(neighbor_list);

}

/*
 * Updates the neighborlist for a particle
 * This is used for all particles that have not moved significantly during an iteration
 *
 * arguments:
 *      particles:      The spatially decomposed particle information
 *      neighbor_list:  The particle's neighborlist
 *      x,y,z:          The indecies for the block which the particle belongs to
 *      array_index:    The index of the particle in the list
 *      moved_list:     The array of lists of which particles have moved
 *      counts:         The number of particles in each of the moved lists
 */
void update_list(Particle* particles[XDiv][YDiv][ZDiv], Neighbor* neighbor_list, int x, int y, int z, int array_index, Particle** moved_list[XDiv][YDiv][ZDiv], int counts[XDiv][YDiv][ZDiv]){

    int j, k;
    Particle *p = &particles[x][y][z][array_index];
    
    for (j = -1; j <= 1; ++j){ //y offset
        for (k = -1; k <= 1; ++k) { //z offset
            check_moved_particles(particles, x+1, y+j, z+k, neighbor_list, *p, moved_list, counts);
        }
    }
    
    check_moved_particles(particles, x, y+1, z+1, neighbor_list, *p,  moved_list, counts);
    check_moved_particles(particles, x, y+1, z  , neighbor_list, *p, moved_list, counts);
    check_moved_particles(particles, x, y+1, z-1, neighbor_list, *p, moved_list, counts);
    check_moved_particles(particles, x, y  , z+1, neighbor_list, *p, moved_list, counts);
    check_moved_particles_center(particles, x, y, z, neighbor_list,  p, moved_list, counts);
    
    sort_list(neighbor_list);
}


/*
 * Builds the entire set of neighborlists for all particles during the first iteration
 *
 * arguments:
 *      particles:  The spatially decomposed particle information
 *      neighbors:  The array of all neighbor lists
 *      my_rank:    The rank of the node
 */
void build_neighbor_lists_initial(Particle* particles[XDiv][YDiv][ZDiv], Neighbor** neighbors, int my_rank){
    int x, y, z, i;
    Particle p;
    
    int tag = 0;
    int rank = my_rank;
    for (x = node_boundries[rank][0]; x <= node_boundries[rank][1]; ++x){
        for (y = node_boundries[rank][2]; y <= node_boundries[rank][3]; ++y){
            for (z = node_boundries[rank][4]; z <= node_boundries[rank][5]; ++z){
                if(in_bounds(x, XDiv) && in_bounds(y, YDiv) && in_bounds(z, ZDiv)){
#pragma omp for private(p)
                    for (i = 0; i < block_size; ++i){
                        p = particles[x][y][z][i];
                        if (p.index != -1) {
                            build_list_full(particles, neighbors[tag + i], x, y, z, i);
                        }
                    }
                    tag = tag + block_size;
                }
            }
        }
    }
}

/*
 * Calculates which particles have moved significantly between iterations, and writes this information out to
 * moved lists
 *
 * arguments:
 *      particles:  The spatially decomposed particle information
 *      moved_lists (OUT):  The array of all particles that have moved
 *      moved_list_counts(OUT): The number of particles in each of the moved lists
 *      my_rank:    The rank of the node
 */
void create_moved_lists(Particle* particles[XDiv][YDiv][ZDiv], Particle** moved_lists[XDiv][YDiv][ZDiv], int moved_list_counts[XDiv][YDiv][ZDiv], int my_rank){
    int x, y, z, i, j;
    Particle* p;
    int index;
    
#pragma omp single
    {
        for (x = 0; x < XDiv; ++x){
            for (y = 0; y < YDiv; ++y) {
                for (z = 0; z < ZDiv; ++z){
                    moved_list_counts[x][y][z] = 0;
                    for (i = 0; i < moved_list_size; ++i) {
                        moved_lists[x][y][z][i] = NULL;
                    }
                }
            }
        }
    }

    for (x = node_boundries[my_rank][0]; x <= node_boundries[my_rank][1]+1; ++x) {
        for (y = node_boundries[my_rank][2]-1; y <= node_boundries[my_rank][3]+1; ++y) {
            for (z = node_boundries[my_rank][4]-1; z <= node_boundries[my_rank][5]+1; ++z) {
                int x_new = (x+XDiv)% XDiv;
                int y_new = (y+YDiv)% YDiv;
                int z_new = (z+ZDiv)% ZDiv;
                if (moved_list_counts[x_new][y_new][z_new] == 0){
                    Particle* p_list = particles[x_new][y_new][z_new];    
#pragma omp for private(p)
                    for (i = 0; i < block_size; ++i){
                        index = 0;
                        p = &p_list[i];
                        if (p->index != -1) {
                            if (p->moved == 1) {
#pragma omp critical
                                {
                                int open_index = moved_list_counts[x_new][y_new][z_new]; 
                                moved_lists[x_new][y_new][z_new][open_index] = p;
                                moved_list_counts[x_new][y_new][z_new]++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


/* Clears the neighborlist */
void clear_list(Neighbor* list){
    int i, j;
    for (i = 0; i < num_neighbors; ++i) {
        list[i].index = -1;
    }
}

/* Removes all the particles that are flagged as moved in the list */
void remove_moved_particles(Neighbor* list){
    int i;
    
    int last;
    for (i = num_neighbors - 1; i >= 0; --i){
        if (list[i].index != -1){
            last = i;
            break;
        }
    }
    for (i = 0; i < num_neighbors; ++i) {
        if (list[i].index != -1 &&                                    
         (list[i].ptr->index != list[i].index ||
         list[i].ptr->moved)){
            if (i != last){
                list[i] = list[last];
                list[last].index = -1;
            }
            else {
                list[i].index = -1;
            }
        }
    }
}


/*
 * Updates all the neighborlists during subsequent iterations
 *
 * arguments:
 *      particles:  The spatially decomposed particle information
 *      neighbors:  The array of all neighbor lists
 *      moved_lists:The list of particles that have moved between iterations
 *      counts:     The number of particles in each moved list
 *      my_rank:    The rank of the node
 */
void update_neighbor_lists(Particle* particles[XDiv][YDiv][ZDiv], Neighbor** neighbors, Particle** moved_list[XDiv][YDiv][ZDiv], int counts[XDiv][YDiv][ZDiv], int my_rank){
    //Iterate through all the particles in all the blocks
    //If a particles has moved, then do a full rebuild of its neighbor list
    //Otherwise, look at the particles that have moved, and update the list
    int x, y, z, i;
    int rank = my_rank - 1;
    int tag = 0;
    Particle p;
    for (x = node_boundries[rank][0]; x <= node_boundries[rank][1]; ++x) {
        for (y = node_boundries[rank][2]; y <= node_boundries[rank][3]; ++y) {
            for (z = node_boundries[rank][4]; z <= node_boundries[rank][5]; ++z) {
                
                #pragma omp for private(p)
                for (i = 0; i < block_size; ++i){
                    p = particles[x][y][z][i];
                    if (p.index != -1) {
                        if (p.moved){
                            clear_list(neighbors[tag+i]);
                            build_list_full(particles, neighbors[tag+i], x, y, z, i);
                        }
                        else {
                            remove_moved_particles(neighbors[tag+i]);
                            update_list(particles, neighbors[tag+i], x, y, z, i, moved_list, counts);
                        }
                    }
                }
                tag = tag + block_size;
            }
        }
    }
}

/*
 *  Receives all the particle information from the master node
 */
void receive_particles(Particle* particles[XDiv][YDiv][ZDiv], int my_rank){
    // Must be single because must be received sequentially
    #pragma omp single
    {
    int x, y, z, i; int indexes;
    int XNew, YNew, ZNew;
    for (x = node_boundries[my_rank][0]; x <= node_boundries[my_rank][1]+1; ++x) {
        for (y = node_boundries[my_rank][2]-1; y <= node_boundries[my_rank][3]+1; ++y) {
            for (z = node_boundries[my_rank][4]-1; z <= node_boundries[my_rank][5]+1; ++z) {
                
                XNew = (x+XDiv)%XDiv;
                YNew = (y+YDiv)%YDiv;
                ZNew = (z+ZDiv)%ZDiv;
                
                indexes = XNew+100*YNew+10000*ZNew;
                MPI_Recv(particles[XNew][YNew][ZNew], block_size, MPI_Particle, 0, indexes, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
    }
}

/* Compresses all the neighbor lists into a single message and sends to the master node
 *
 * arguments:
 *      neighbors: The set of neighbor lists that are being sent
 *      particles: The spatially decomposed particle information
 *      my_rank:   Rank of the node
 *      nvlists:   The arrays which mark where each list ends
 */
void send_neighborlists(Neighbor** neighbors, Particle* particles[XDiv][YDiv][ZDiv], int my_rank, int* nvlists){
    
#pragma omp single
    {
    int x, y, z, i;
    int j = 0;
    int sum[NUM_LISTS];

    // Allocate the message array, which is used to send the lists
    int size = num_neighbors * (n / num_work_nodes);
    int * message = malloc(sizeof(int) * size);
    int k = 0;
    int index = 0;
    int temp_index;
    
        
    // Iterate through all the particle arrays, and for each particle in the blocks
    // copy its neighbor list to the message array
    for (x = node_boundries[my_rank][0]; x <= node_boundries[my_rank][1]; ++x) {
        for (y = node_boundries[my_rank][2]; y <= node_boundries[my_rank][3]; ++y) {
            for (z = node_boundries[my_rank][4]; z <= node_boundries[my_rank][5]; ++z) {
                for (i = 0; i < block_size; ++i) {
                    if (particles[x][y][z][i].index != -1) {
                        temp_index = particles[x][y][z][i].index;
                        message[index] = -1; index++;
                        message[index] = temp_index; index++;
                        
                        //The sum array is used to calculate the nvlst array, which is
                        // where each list ends
                        for (j = 0; j < NUM_LISTS; ++j){
                            sum[j] = 0;
                        }
                        for (j = 0; j < num_neighbors; ++j) {
                            if(neighbors[k][j].index != -1) {
                                sum[neighbors[k][j].list - 1] += 1;
                                message[index] = neighbors[k][j].index; index++;
                            }
                        }
                        
                        for (j = 0; j < NUM_LISTS; ++j){
                            if (j != 0) sum[j] = sum[j-1] + sum[j];
                            nvlists[j*n + pIndex(temp_index)] = sum[j];
                        }
                    }
                    k++;
                }
            }
        }
    }
    message[index] = -2; index++;
    
    MPI_Send(message, index , MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(nvlists, NUM_LISTS * n, MPI_INT, 0, 0, MPI_COMM_WORLD);
    free(message);
    }
    

    
}


/*
 * Recieves an instruction from the master node, and write it out to value
 */
void get_instructions(int my_rank, int* value){
#pragma omp single
    {
    MPI_Recv(value, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}


/*
 * Allocates and initializes the data structures used in a neighbor list build by the work nodes
 * 
 * arguments:
 *      num_particles:  number of particles in the system
 *      my_rank:        rank of the node
 *      comm_sz:        number of nodes
 *      neighbors:      neighbor list array
 *      particles:      the arrays holding the spatially decomposed particle arrays
 *      moved_list:     the arrays holding the lists of particles that have moved singificantly
 *      moved_list_counts: the number of particles in each moved list
 */
void initialize_data_structures_work(int num_particles,  int my_rank, int comm_sz, Neighbor*** neighbors, Particle* particles[XDiv][YDiv][ZDiv], Particle**moved_list[XDiv][YDiv][ZDiv], int moved_list_counts[XDiv][YDiv][ZDiv]){
    int x, y, z, i, j;
    int xwidth, ywidth, zwidth;
    
#pragma omp single
    {
        block_size = max(20, 2*n/(XDiv * YDiv * ZDiv));
        moved_list_size = .2 * block_size;
        determine_boundries(comm_sz, my_rank);
    }

    if(node_boundries[my_rank-1][0] != -1){   
#pragma omp single
        {
            xwidth = (node_boundries[my_rank-1][1] - node_boundries[my_rank-1][0]) +1; //The + 1 is because the boundry limits are inclusive
            ywidth = (node_boundries[my_rank-1][3] - node_boundries[my_rank-1][2]) +1;
            zwidth = (node_boundries[my_rank-1][5] - node_boundries[my_rank-1][4]) +1;
            maxnlist = xwidth * ywidth * zwidth * block_size;
            *neighbors = malloc(maxnlist  * sizeof(Neighbor*));
        }

#pragma omp for
        for (i = 0; i < maxnlist; ++i){
            (*neighbors)[i] = malloc(sizeof(Neighbor) * num_neighbors); 
            for (j = 0; j < num_neighbors; ++j){
                (*neighbors)[i][j].index = -1;
            }
        }

#pragma omp for        
        for (x = 0; x < XDiv; ++x){
            for (y = 0; y < YDiv; ++y){
                for (z = 0; z < ZDiv; ++z){
                   particles[x][y][z] = malloc(sizeof(Particle) * block_size);
                   moved_list[x][y][z] = malloc(sizeof(Particle*) * moved_list_size);
                }
            }
        }

#pragma omp for            
        for (x = 0; x < XDiv; ++x) {
            for (y = 0; y < YDiv; ++y) {
                for (z = 0; z < ZDiv; ++z) {
                    moved_list_counts[x][y][z] = 0;
                }
            }
        }

#pragma omp for
        for (i = 0; i < NUM_LISTS; ++i){
            THRESHOLDS_2[i] = THRESHOLDS[i]*THRESHOLDS[i];
        }
        
        nlists = malloc(sizeof(int) * NUM_LISTS * n);
    }
}

/*
 * The initial neighbor list build for the work node
 *
 * arguments:
 *      particles:  The array that holds the spatially decomposed particle information
 *      my_rank:    The rank of the node
 *      neighbors:  The array that will hold the neighborlists
 */

void initial_build_work(Particle* particles[XDiv][YDiv][ZDiv], int my_rank, Neighbor** neighbors){
    receive_particles(particles, my_rank-1);
    build_neighbor_lists_initial(particles, neighbors, my_rank - 1);
    send_neighborlists(neighbors, particles, my_rank - 1, nlists);
}

/*
 * Every list build following the initial one
 *
 * arguments:
 *     particles: The array that holds the spatially decomposed particle information
 *     my_rank:   The rank of the node
 *     moved_list_counts:   The number of particles in each moved list
 *     moved_lists: The array that will hold the list of particles that have moved between iterations
 *     neighbors:   The array that contains the neighborlists from the prior iteration
 */
void subsequent_build_work(Particle* particles[XDiv][YDiv][ZDiv], int my_rank, int moved_list_counts[XDiv][YDiv][ZDiv], Particle** moved_list[XDiv][YDiv][ZDiv], Neighbor** neighbors){
    receive_particles(particles, my_rank-1);
    create_moved_lists(particles, moved_list, moved_list_counts, my_rank - 1);
    update_neighbor_lists(particles, neighbors, moved_list, moved_list_counts, my_rank);
    send_neighborlists(neighbors, particles, my_rank - 1, nlists);
}


/*
 * Frees memory used by the work nodes at the end
 */
void clean_up_work(Particle* particles[XDiv][YDiv][ZDiv], Particle** moved_list[XDiv][YDiv][ZDiv], Neighbor** neighbors){
    int x, y, z, i;
#pragma omp for
    for (x = 0; x < XDiv; ++x){
        for (y = 0; y < YDiv; ++y){
            for (z = 0; z < ZDiv; ++z){
                free(particles[x][y][z]);
                free(moved_list[x][y][z]);
            }
        }
    }
    
#pragma omp for
    for (i = 0; i < maxnlist; ++i){
        free(neighbors[i]);
    }
    
#pragma omp single
    {
        free(neighbors);
        free(nlists);
    }
    
}

/*
 * Recieves various parameters that remain constant throughout the simulations
 * needed by the work node to build the list
 */
void receive_extra_information(){
#pragma omp single 
    {
    int int_info[6];
    double double_info[3];
        
    MPI_Recv(int_info, 6, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(double_info, 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(THRESHOLDS, 10, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
    XDiv = int_info[0];
    YDiv = int_info[1];
    ZDiv = int_info[2];
    n    = int_info[3];
    num_neighbors   = int_info[4];
    NUM_LISTS = int_info[5];
    
    DIV_WIDTH_X = double_info[0];
    DIV_WIDTH_Y = double_info[1];
    DIV_WIDTH_Z = double_info[2];
    
    WIDTH_X = XDiv * DIV_WIDTH_X;
    WIDTH_Y = YDiv * DIV_WIDTH_Y;
    WIDTH_Z = ZDiv * DIV_WIDTH_Z;
        
    WIDTH_X2 = WIDTH_X/2;
    WIDTH_Y2 = WIDTH_Y/2;
    WIDTH_Z2 = WIDTH_Z/2;
    }
}




