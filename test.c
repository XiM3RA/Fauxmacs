#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define DIM 2
typedef double real;
#define sqr(x) ((x)*(x))

typedef struct {
    real m;
    real x[DIM];
    real v[DIM];
    real F[DIM];
    real F_old[DIM];
} Particle;

void timeIntegration_basis(real t, real delta_t, real t_end,
                           Particle *p, int N);
void compoutStatistic_basis(Particle *p, int N, real t);
void compX_basis(Particle *p, int N, real delta_t);
void compV_basis(Particle *p, int N, real delta_t);
void updateX(Particle *p, real delta_t);
void updateV(Particle *p, real delta_t);
void compF_basis(Particle *p, int N);
void outputResults_basis(Particle *p, int N, real t);
void force(Particle *i, Particle *j);
void inputParameters_basis(real delta_t, real t_end, int N);
void getRunParameters(double* delta_t, double* t_end, int *N);

void force(Particle *i, Particle *j) {
    real r = 0;
    for (int d = 0; d < DIM; d++)
        r += sqr(j->x[d] - i->x[d]);
    real f = i->m * j->m / (sqrt(r) * r);
    for (int d = 0; d < DIM; d++)
        i->F[d] += f * (j->x[d] - i->x[d]);
}

void outputResults_basis(Particle *p, int N, real t) {
    printf("%lf\n", t);
    int i = 0;
    while(i < N) {
        printf("%lf %lf %lf %lf\n", p[i].x[0], p[i].x[1], p[i].v[0], p[i].v[1]);
        ++i;
    }
}

void timeIntegration_basis(real t, real delta_t, real t_end,
                           Particle *p, int N) {
    compF_basis(p, N);
    while (t < t_end) {
        t += delta_t;
        compX_basis(p, N, delta_t);
        compF_basis(p, N);
        compV_basis(p, N, delta_t);
        compoutStatistic_basis(p, N, t);
        // Need an output function
        outputResults_basis(p, N, t);
    }
}

void compoutStatistic_basis(Particle *p, int N, real t) {
    real e = 0;
    for (int i=0; i < N; i++) {
        real v = 0;
        for (int d = 0; d < DIM; d++)
            v += sqr(p[i].v[d]);
        e += 0.5 * p[i].m * v;
    }
}

void compX_basis(Particle *p, int N, real delta_t) {
    for (int i = 0; i < N; i++)
        updateX(&p[i], delta_t);
}

void compV_basis(Particle *p, int N, real delta_t) {
    for (int i = 0; i < N; i++)
        updateV(&p[i], delta_t);
}

void compF_basis(Particle *p, int N) {
    for (int i = 0; i < N; i++)
        for (int d = 0; d < DIM; d++)
            p[i].F[d] = 0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            if (i != j) force(&p[i], &p[j]);
}

void updateX(Particle *p, real delta_t) {
    real a = delta_t * 0.5 / p->m;
    for (int d = 0; d < DIM; d++) {
        p->x[d] += delta_t * (p->v[d] + a * p->F[d]);
        p->F_old[d] = p->F[d];
    }
}

void updateV(Particle *p, real delta_t) {
    real a = delta_t * 0.5 / p->m;
    for (int d = 0; d < DIM; d++)
        p->v[d] += a * (p->F[d] + p->F_old[d]);
}

void getRunParameters(double* delta_t, double* t_end, int *N) {
    printf("Enter timestep, trajectory length, number of particles (separated by spaces): ");
    scanf("%lf %lf %o", delta_t, t_end, N);
}


int main () {
    int N;
    double delta_t, t_end;
    char *filename = "run-params.txt";
    FILE *fp = fopen(filename, "r");

    if (fp == NULL)
    {
        printf("Error: could not open file %s", filename);
        return 1;
    }

    fscanf (fp, "%d %lf %lf", &N, &delta_t, &t_end);
    Particle *p = (Particle*)malloc(N * sizeof(*p));

    int i = 0;
    while (!feof (fp)) {
        fscanf (fp, "%lf %lf %lf %lf %lf", &p[i].m, &p[i].x[0], &p[i].x[1],
                                                &p[i].v[0], &p[i].v[1]);
        ++i;   
    }

    timeIntegration_basis(0, delta_t, t_end, p, N);
    free(p);
    return 0;
}