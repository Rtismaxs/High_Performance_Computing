#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>

// Function that executes the branch-and-bound algorithm
// INPUT: total number of cities (int), distance matrix (array), current path (array), cities visited (array), bound (int), minimal cost (int), depth (int), path with the minimal cost (array)
// OUTPUT: None
void wsp(int num_cities, int distance_matrix[num_cities+1][num_cities+1], int path[num_cities], int visited[num_cities], int bound, int *min_tour, int depth, int final_path[num_cities])
{
    
    // When the bound is superior to the minimal cost, we stop the current path
    if (bound > *min_tour){
        return;
    }
    
    // When the depth is equal to the number of cities or when we visited all the cities
    if (depth == num_cities){
        // If the current bound is inferior to the minimal cost
        if (bound < *min_tour){
            // Change the value of the minimal cost
            *min_tour = bound;
            // Store the specific path
            for (int i = 0; i < num_cities; i++) {
                final_path[i] = path[i];
            }
        }
        return;
    }
    
    // Go through all the cities
    for (int city = 1; city < num_cities+1 ; city++) {
        // If we have not already visited the city
        if (!visited[city-1]){
            
            // Add the city to the current path and update the visited array
            path[depth]=city;
            visited[city-1] = 1;
            
            // Variable to store the value of the current bound 
            int temp = bound;
            
            // Update the bound by adding the distance between the previous city and the current city
            bound += distance_matrix[path[depth-1]][city];
            
            // Recursivity and change in the depth at which the path is continued
            wsp(num_cities, distance_matrix, path, visited, bound, min_tour, depth + 1, final_path);
            
            // Recover the previous value of the bound before the adding up and the recursivity
            bound = temp;
            // Update the visited array for backtracking
            visited[city-1]=0;
        }
    }
    return;
}

// Function that swaps the memory from 2 pointers
// INPUT: int, int
// OUTPUT: None
void swap(int *x, int *y) {
    int temp = *x;
    *x = *y;
    *y = temp;
}

// Global variable used as a counter for the calculation of the permutations
int counter = 0;

// Function that generates a specific permutation with a specific length
// INPUT: total number of cities (int), permutation generator (array), permutation storage (array), permutation from which integer (int), permutation until which integer (int), total length of the permutation (int), first integer of the permutation (int), counter (int)
// OUTPUT: None
void permute(int num_cities, int *lst, int permute_path[num_cities], int start, int end, int rank, int length, int random_number) {
    
    // When the index of the permutation is equal to the specified length
    if (start == length) {
    
        // Keep only the permutation that begins with the specified first integer
        if(lst[0] == random_number) {
        
            // Store only the specific permutation that corresponds to the rank specified
            if (rank == counter){
            
                // Save the appropriate permutation 
                for (int i = 0; i < length; i++) {
                    permute_path[i] = lst[i];
                }
                
                // Add 0 to have a length of the pemutation equal to the total number of cities
                for (int i = length; i < num_cities; i++) {
                    permute_path[i] = 0;
                }
                // Add 1 to the counter to continue with the next possible permutation
                counter +=1;
                
            } else {
                
                // Add 1 to the counter to continue with the next possible permutation
                counter +=1;
            }
        }
    } else {
        // Calculate permutations from a certain index until another specified index 
        for (int i = start; i <= end; i++) {
            // Operate swaping
            swap((lst + start), (lst + i));
            // Recursivity and change in the start at which the permutation is continued
            permute(num_cities, lst, permute_path, start + 1, end, rank, length, random_number);
            // Operate swaping after finishing permutations from a certain start
            swap((lst + start), (lst + i));
        }
    }
}

// Main function
void main(int argc, char *argv[])
{
    
    // MPI Initialisation
    MPI_Init(&argc, &argv);
    
    // Variables for the number of process and the rank of the specific process
    int npes, myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
	  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    // Variables to measure the time
    double start, end;
    
    // Start the timer
    start = MPI_Wtime();

    // Recover the input file with the -i command
  	char* input_file = NULL;
  	for (int i = 0; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0) {
   		      input_file = argv[i+1];
  	   	}
	  }
     
    // Case where the input file is not specified
  	if (input_file == NULL) {
  	    fprintf(stderr, "Error: input file not specified.\n");
   		  return;
	  }
    
    // We open the input file
  	FILE* file = fopen(input_file, "r");
    
     // Case where we can't open the input file
	  if (file == NULL) {
        fprintf(stderr, "Error: unable to open input file.\n");
        return;
	  }
    
    // Variable for the total number of cities
    int num_cities;
    // Read the total number of cities from the first line of the input file
	  fscanf(file, "%d", &num_cities);
    
    // Variable for the distance matrix
    int distance_matrix[num_cities+1][num_cities+1];
    // Read the distance matrix from the input file
	  for (int i = 1; i < num_cities+1; i++) {
	      for (int j = 1; j < i; j++) {
		        fscanf(file, "%d", &distance_matrix[i][j]);
            // We want a full matrix so distance_matrix[j][i] = distance_matrix[i][j]
            distance_matrix[j][i] = distance_matrix[i][j];
        }
        // We add 0 on the first diagonal
        distance_matrix[i][i] = 0;
    }
    
    // We close the input file
    fclose(file);
    
    // Display the total number of cities and the distance matrix with the process 0
    if (myrank == 0){
        printf("For %d cities with the respective distance matrix:", num_cities);
        for (int i = 1; i < num_cities+1; i++) {
                for (int j = 1; j < i; j++) {
                    printf("%d ",distance_matrix[i][j]);
            }
            printf("\n");
        }
        printf("\n");
    }
    
    // Initialisation of the minimal cost with a high value
    int min_tour = INT_MAX;

    // Variable for the bound of the path
    int bound = 0;
    // Variable for the depth of the path
    int depth;
    
    // Variable to store the path with the minimal cost
    int final_path[num_cities];
    // Variable to store the permutations
    int permute_path[num_cities];
    // Variable to store the path we take
    int path[num_cities];
    // Variable to store cities we already visited
    int visited[num_cities];
    
    // Variable for calculating the quotient when there more or less process than the total number of cities
    int quotient;
    // Variable for calculating the remainder when there more or less process than the total number of cities
    int remainder;

    // Variable to store the minimal cost and rank for each process
    int root[2];
    // Variable to store the minimal cost and rank where it is stored
    int result[2];

    // Generate a random number between 1 and num_cities
    srand(time(NULL));
    int random_number = rand() % num_cities + 1;

    // Variable to store the number of cities - 1
    int num_curr_cities = num_cities-1;
    // Variable to store the number of nodes for a specific depth
    int f = num_curr_cities;

    // We start at a depth of 2
    depth = 2;
    
    // Calculate the quotient between the number of nodes on the first depth and the number of processors   
    quotient = (f)/(npes);
    
    // While loop to find the best depth where the quotient is higher than the number of processors
    // The idea is to refer to the case where we have fewer processors than nodes
    while((num_curr_cities > 1) && (quotient < npes)){
        depth += 1;
        num_curr_cities -= 1;
        f = f * num_curr_cities;
        quotient = (f)/(npes);    
    }
    
    // Calculate the remainder between number of nodes on a certain depth and the number of processors
    remainder = (f)%(npes);

    // Case where there are the same number or more processors than total number of possible paths
    if (npes >= f){
        
        // Create a new communicator that contains only the first processes until the total number of paths
        MPI_Comm new_comm;
        
        if (myrank < f) {
            MPI_Comm_split(MPI_COMM_WORLD, 0, myrank, &new_comm);
        } else {
            MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, myrank, &new_comm);
        }
        
        // Going trough all the others process until the total number of paths
        for (int i=0; i < f; i++){
            if (myrank == i){
                
                // Initialisation of these two arrays: path and visited
                for (int j = 0; j < num_cities; j++) {
                    path[j]=j+1;
                    visited[j]=0;
                }
                
                // Calculate all the permutations possible for a depth of num_cities starting with the random number generated
                // Store the permutation n°i into the permute_path array
                permute(num_cities, path, final_path, 0, num_cities - 1, i, num_cities, random_number);
                                    
                // Reinitialise the counter to 0
                counter = 0;
                
                // Initialise the bound to 0
                bound = 0;

                // Calculate the initial bound based on the permutation we generated
                for (int i = 0; i < num_cities-1; i++) {
                    bound += distance_matrix[final_path[i]][final_path[i+1]];
                }
                
                min_tour = bound;
                
                // Associate the minimal cost and rank
                root[0] = min_tour;
                root[1] = myrank;
                
                // Apply the MPI_Allreduce method to find the minimal cost trough all the processes and to send this value to all the processes
                MPI_Allreduce(MPI_IN_PLACE, &min_tour, 1, MPI_INT, MPI_MIN, new_comm);

                // Apply the MPI_Allreduce method to find the minimal cost and its location trough all the processes and to send these values to all the processes
                MPI_Allreduce(root, result, 1, MPI_2INT, MPI_MINLOC, new_comm);
                
                // Broadcast of the path with the minimal cost                  
                MPI_Bcast(final_path, num_cities, MPI_INT, result[1], new_comm);                                                     
            }
        }
    }

    // Case where there are less processors than nodes on a certain depth
    else if (npes < f){
        
        // Case where there is only 1 worker
        if (quotient == f){
            
            // Going trough all the others process
            for (int i=0; i < npes; i++){
                if (myrank == i){
                    
                    // Initialisation of these two arrays: path and visited
                    for (int j = 0; j < num_cities; j++) {
                        path[j]=0;
                        visited[j]=0;
                    }
                    
                    // Fix the first element of the path with this city and update the visited array
                    path[0] = random_number;
                    visited[path[0]-1] = 1;
                    
                    // In this case, we are in depth of 1
                    depth = 1;
                    
                    // Initialise the bound to 0
                    bound = 0;
                    
                    // Run the branch-and-bound algorithm
                    wsp(num_cities, distance_matrix, path, visited, bound, &min_tour, depth, final_path);      
                }
            }
            
        // Case where the quotient is not equal to the number of nodes on a certain depth
        } else if (quotient != f){
    
            int c = 0;
            // Going trough all the others process as many times as the value of the quotient calculated previously
            for (int q=0; q < quotient; q++){
                for (int i=0; i < npes; i++){
                    if (myrank == i){
                    
                        // Initialisation of these two arrays: path and visited
                        for (int j = 0; j < num_cities; j++) {
                            path[j]=j+1;
                            visited[j]=0;
                        }
                    
                        // Calculate all the permutations possible for a depth of num_cities starting with the random number generated
                        // Store the permutation n°i into the permute_path array
                        permute(num_cities, path, permute_path, 0, num_cities - 1, i+c, depth, random_number);
                                        
                        // Reinitialise the counter to 0
                        counter = 0;
                            
                        // Initialise the bound to 0
                        bound = 0;
                    
                        // Calculate the initial bound based on the permutation we generated
                        for (int i = 0; i < depth-1; i++) {
                            bound += distance_matrix[permute_path[i]][permute_path[i+1]];
                        }
                        
                        // Update the path and visited array based on the permutation we generated
                        for (int j = 0; j < num_cities; j++) {
                            path[j]=permute_path[j];
                            visited[path[j]-1]=1;
                        }
                       
                        // Run the branch-and-bound algorithm
                        wsp(num_cities, distance_matrix, path, visited, bound, &min_tour, depth, final_path);
                        
                        // Associate the minimal cost and rank
                        root[0] = min_tour;
                        root[1] = myrank;
                        
                        // Apply the MPI_Allreduce method to find the minimal cost trough all the processes and to send this value to all the processes
                        MPI_Allreduce(MPI_IN_PLACE, &min_tour, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
                        
                        // Apply the MPI_Allreduce method to find the minimal cost and its location trough all the processes and to send these values to all the processes
                        MPI_Allreduce(root, result, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);
                        
                    }
                }
                // Add npes to the counter that keeps us informed about the number of times we have gone through the quotient loop
                c += npes;
                
                // Broadcast of the path with the minimal cost 
                MPI_Bcast(final_path, num_cities, MPI_INT, result[1], MPI_COMM_WORLD);        
            }
    
            // Case where the remainder is not null, we are going trough all processes as many times as the value of the remainder calculated previously
            if (remainder != 0){
                
                // Create a new communicator that contains only the first remainder processes
                MPI_Comm new_comm;
                
                if (myrank < remainder) {
                    MPI_Comm_split(MPI_COMM_WORLD, 0, myrank, &new_comm);
                } else {
                    MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, myrank, &new_comm);
                }
            
                // Going trough all the others process as many times as the value of the remainder calculated previously
                for (int r=0; r < remainder; r++){
                
                    if (myrank == r){
                        
                        // Initialisation of these two arrays: path and visited
                        for (int j = 0; j < num_cities; j++) {
                            path[j]=j+1;
                            visited[j]=0;
                        }
                        
                        // Calculate all the permutations possible for a depth of num_cities starting with the random number generated
                        // Store the permutation n°i into the permute_path array
                        permute(num_cities, path, permute_path, 0, num_cities - 1, f-myrank-1, depth, random_number); 
                                            
                        // Reinitialise the counter to 0
                        counter = 0;
                            
                        // Initialise the bound to 0
                        bound = 0;
                        
                        // Calculate the initial bound based on the permutation we generated
                        for (int i = 0; i < depth-1; i++) {
                            bound += distance_matrix[permute_path[i]][permute_path[i+1]];
                        }
                        
                        // Update the path and visited array based on the permutation we generated
                        for (int j = 0; j < num_cities; j++) {
                            path[j]=permute_path[j];
                            visited[path[j]-1]=1;
                        }
                        
                        // Run the branch-and-bound algorithm
                        wsp(num_cities, distance_matrix, path, visited, bound, &min_tour, depth, final_path);
                        
                        // Associate the minimal cost and rank
                        root[0] = min_tour;
                        root[1] = myrank;
                        
                        // Apply the MPI_Allreduce method to find the minimal cost trough all the processes and to send this value to all the processes
                        MPI_Allreduce(MPI_IN_PLACE, &min_tour, 1, MPI_INT, MPI_MIN, new_comm);
                        
                        // Apply the MPI_Allreduce method to find the minimal cost and its location trough all the processes and to send these values to all the processes
                        MPI_Allreduce(root, result, 1, MPI_2INT, MPI_MINLOC, new_comm);
                        
                        // Broadcast of the path with the minimal cost
                        MPI_Bcast(final_path, num_cities, MPI_INT, result[1], new_comm);
                        
                    }
                }
            }
        }               
    }
                            
    // Display the results obtained from the parallelised branch-and-bound algorithm
    if (myrank == 0){
        
        // End the timer
        end = MPI_Wtime();
        
        printf("\n");
        // Display the random city we generated and started with
        printf("We will start with the city %d\n", random_number);
        
        // Display the path with the minimal cost
        printf("We obtained the following path as a solution of the WSP: ");
        for (int i = 0; i < num_cities; i++) {
            printf("%d ", final_path[i]);
        }
        printf("\n");
        
        // Display the minimal cost calculated with the branch-and-bound algorithm
        printf("The minimal bound calculated is: %d \n", min_tour);
        
        // Calculate and display the runtime in seconds
        printf("It takes %f seconds to calculate the Branch and Bound algorithm\n", end-start);
    }            
    
    // Finalise MPI
    MPI_Finalize();
    
    return;
}