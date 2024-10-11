#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

/*
Populates the given vector/matrix */
void populate_x(int *x);    //klar
void populate_p(int **p);   //klar

// Given propensity functions
void prob(int *x, double *w); // klar

// Stochastic Simulation algorithms
// int SSA(int **p, int **local_results, int iter, int thread_id, int *x, double *w);
int SSA(int **p, int **results, int iter);
void update_state_vector(int *x, int **p, int r);

// General helper functions
double sum_w(double *w, int len);
int find_r(double *w, double a0u2);
void update_state_vector(int *x, int **p, int r);

// function to create output
void create_output(const char *filename, int **matrix, int cols);
