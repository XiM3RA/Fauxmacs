#ifndef UTILITIES
#define UTILITIES

#define DIM 2
typedef double real;

typdef struct {
    real m;
    real x[DIM];
    real v[DIM];
    real F[DIM];
} Particle;

#endif