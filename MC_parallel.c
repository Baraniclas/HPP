#include "MC.h"

// Take N as argument (num of experiments)

int main(int argc, char **argv) {

    if (argc < 2) {
        printf("Expected 1 argument, number of simulations, N \n");
        exit(1);
    }

    int N = atoi(argv[1]);
    if (N <= 0) {
        printf("N must be a positive integer. \n");
        exit(1);
    }

    // allocate memory for matrix P
    int **p = (int **)malloc(15*sizeof(int*));

    if (p == NULL) {
        printf("Memory allocation for matrix P failed. \n");
        exit(1);
    }

    for (int ind = 0; ind < 15; ind ++) {
        p[ind] = (int *)calloc(7, sizeof(int));
    }



    int ** results = (int **)malloc(7*sizeof(int*));
    if (results == NULL) {
       printf("Memory allocation for results matrix failed. \n");
       exit(1);
   }      
    for (int ind = 0; ind < 7; ind++) {
        results[ind] = (int *)malloc(N*sizeof(int));
    }
    
    // Populate matrix P gives previously known values.
    populate_p(p);

    // begin parallelization with dynamic scheduling.
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < N; i ++) {
            SSA(p, results, i);
        }
    }

    // generates the output for the simulations within a CSV file.
    create_output("ResultsMatrix.csv", results, N);

    for (int ind = 0; ind < 7; ind++) {
       free(results[ind]);
   }
   free(results);

   for (int ind = 0; ind < 15; ind++) {
       free(p[ind]);
   }
   free(p);

    return 0;
}

int SSA(int **p, int **results, int iter) {
    int thread_id = omp_get_thread_num();

    unsigned int seed = 5*(iter+1)*thread_id+10000;

    double end_time = 100.0, curr_time = 0.0;
    double u1, u2, tau, a0, a0u2;
    int r;

    int *x = (int *)malloc(7*sizeof(int));
    double *w = (double *)malloc(15*sizeof(double));


    if (w == NULL) {
        printf("Memory allocation of vector w failed. \n");
        return 1;
    }
    if (x == NULL) {
        printf("Memory allocation of vector x failed. \n");
        return 1;
    }

    populate_x(x);

    while (curr_time < end_time) {
        prob(x, w);

        a0 = sum_w(w, 15);

        u1 = (double)(rand_r(&seed) + 1.0) / (RAND_MAX + 1.0);
        u2 = (double)(rand_r(&seed) + 1.0) / (RAND_MAX + 1.0);
        
        tau = -log(u1)/a0;
        a0u2 = a0*u2;

        r = find_r(w, a0u2);   

        update_state_vector(x, p, r);

        curr_time += tau;
    }

    // enter results from SSA run into local_results
    for (int j = 0; j < 7; j++){
        results[j][iter] = x[j];
    }

    free(x);
    free(w);

    return 0;
}

void update_state_vector(int *x, int **p, int r) {
    for (int i = 0; i < 7; i++) {
        x[i] = x[i] + p[r][i];
    }
}

int find_r(double *w, double a0u2) {
    double x, y;
    y = 0;
    for (int r = 0; r < 15; r++) {
        x = y;
        y += w[r];
        if (x < a0u2 && a0u2 <= y) {
            return r;
        }
    }
    return 0;
}

double sum_w(double *w, int len) {
    double sum = 0.0;
    for (int i = 0; i < len; i++) {
        sum += w[i];
    }
    return sum;
}

void populate_x(int *x) {    
    x[0] = 900;
    x[1] = 900;
    x[2] = 30;
    x[3] = 330;
    x[4] = 50;
    x[5] = 270;
    x[6] = 20;
}

void populate_p(int **p) {
    p[0][0] = 1;
    p[1][0] = -1;
    p[2][0] = -1;
    p[2][2] = 1;
    p[3][1] = 1;
    p[4][1] = -1;
    p[5][1] = -1;
    p[5][3] = 1;
    p[6][2] = -1;
    p[7][2] = -1;
    p[7][4] = 1;
    p[8][3] = -1;
    p[9][3] = -1;
    p[9][5] = 1;
    p[10][4] = -1;
    p[11][4] = -1;
    p[11][6] = 1;
    p[12][5] = -1;
    p[13][0] = 1;
    p[13][6] = -1;
    p[14][6] = -1;
}

/**
 * Compute propensities for the Malaria model.
 * @param x State vector Should be of length 7!
 * @param w Result vector (propensities). Should be of length 15!
 *
 */
void prob(int *x, double *w) {
	
	const double LAMBDA_H = 20;     // Birth number, humans
	const double LAMBDA_M = 0.5;    // Birth number, mosquitoes
	const double B = 0.075; 	    // Biting rate of mosquitoes
	const double BETA_H = 0.3;  	// Probability that a bite by an infectious mosquito results in transmission of disease to human
	const double BETA_M = 0.5;  	// Probability that a bite results in transmission of parasite to a susceptible mosquito
	const double MU_H = 0.015;      // Human mortality rate
	const double MU_M = 0.02;       // Mosquito mortality rate
	const double DELTA_H = 0.05;    // Disease induced death rate, humans
	const double DELTA_M = 0.15;    // Disease induced death rate, mosquitoes
	const double ALFA_H = 0.6;      // Rate of progression from exposed to infectious state, humans
	const double ALFA_M = 0.6;      // Rate of progression from exposed to infectious state, mosquitoes
	const double R = 0.05;          // Recovery rate, humans
	const double OMEGA = 0.02;      // Loss of immunity rate, humans
	const double NU_H = 0.5;        // Proportion of an antibody produced by human in response to the incidence of infection caused by mosquito
	const double NU_M = 0.15;       // Proportion of an antibody produced by mosquito in response to the incidence of infection caused by human

	w[0] = LAMBDA_H;
	w[1] = MU_H * x[0];
	w[2] = (B * BETA_H * x[0] * x[5]) / (1 + NU_H * x[5]);
	w[3] = LAMBDA_M;
	w[4] = MU_M * x[1];
	w[5] = (B * BETA_M * x[1]*x[4]) / (1 + NU_M * x[4]);
	w[6] = MU_H * x[2];
	w[7] = ALFA_H * x[2];
	w[8] = MU_M * x[3];
	w[9] = ALFA_M * x[3];
	w[10] = (MU_H + DELTA_H) * x[4];
	w[11] = R * x[4];
	w[12] = (MU_M + DELTA_M) * x[5];
	w[13] = OMEGA * x[6];
	w[14] = MU_H * x[6];
}

// this function is based off of code from this forum:
// https://stackoverflow.com/questions/42550882/write-values-from-multidimensional-array-to-csv-file
void create_output(const char *filename, int **matrix, int cols) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error opening file %s\n", filename);
        return;
    }
    for (int i = 0; i < 7; i ++) {
        for (int j = 0; j < cols; j++) {
            fprintf(file, "%d%s", matrix[i][j], (j < cols - 1) ? "," : "");
        }
        fprintf(file, "\n");
    }
    fclose(file);
}