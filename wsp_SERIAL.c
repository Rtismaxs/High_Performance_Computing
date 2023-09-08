#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <sys/time.h>

// Function that executes the branch-and-bound algorithm
void wsp(int num_cities, int distance_matrix[num_cities+1][num_cities+1], int path[num_cities], int visited[num_cities], int bound, int *min_tour, int depth, int final_path[num_cities]);

// Main function
void main(int argc, char *argv[])
{
    // Variables to measure the time
    struct timeval start, end;
    double runtime;

    // Start the timer
    gettimeofday(&start, NULL);
    
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
    
    // We display the total number of cities
    printf("For %d cities with the respective distance matrix:", num_cities);
    
    // We display the distance matrix as a lower triangular matrix
    for (int i = 1; i < num_cities+1; i++) {
	      for (int j = 1; j < i; j++) {
		        printf("%d ",distance_matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    
    // Variable to store the path we take
    int path[num_cities];
    // Variable to store cities we already visited
    int visited[num_cities];
    
    // Initialisation of these two arrays
    for (int i = 0; i < num_cities; i++) {
        path[i]=0;
        visited[i]=0;
    }

    // Generate a random number between 1 and num_cities
    //srand(time(NULL));
    //int random_number = rand() % num_cities + 1;
    int random_number = 1;
    // Display the random city we generated
    printf("We will start with the city %d\n", random_number);

    // Fix the first element of the path with this city and update the visited array
    path[0]=random_number;
    visited[path[0]-1]=1;
    
    // Initialisation of the bound of the path
    int bound = 0;
    
    // Initialisation of the minimal cost with a high value
    int min_tour = INT_MAX;
    
    // Initialisation of the depth we start the branch-and-bound algorithm
    int depth = 1;
    
    // Variable to store the path with the minimal cost
    int final_path[num_cities];
    
    // Run the branch-and-bound algorithm
    wsp(num_cities, distance_matrix, path, visited, bound, &min_tour, depth, final_path);
    
    // End the timer
    gettimeofday(&end, NULL);
    
    // Display the path with the minimal cost
    printf("We obtained the following path as a solution of the WSP: ");
    for (int i = 0; i < num_cities; i++) {
        printf("%d ", final_path[i]);
    }
    printf("\n");
    
    // Display the minimal cost calculated with the branch-and-bound algorithm
    printf("The minimal bound calculated is: %d \n", min_tour);
    
    // Calculate the runtime in seconds
    runtime = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    // Display the runtime
    printf("It takes %f seconds to calculate the Branch and Bound algorithm\n", runtime);
    
    return;
}

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