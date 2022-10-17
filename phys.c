#include <math.h>

#include <stdio.h>

#include "phys.h"

/**
 * Lambda in cgs
*/
double planck(double temp, double lambda) {
    return 2 * h * pow(c, 2.) / (pow(lambda, 5.) * (exp(h * c / (lambda * k * temp)) - 1));
}

double e1(double x) {
    if (x > 0. && x <= 1)
        return -1. * log(x) -
                     0.57721566 +
                     0.99999193 * x -
                     0.24991055 * pow(x, 2.) +
                     0.05519968 * pow(x, 3.) -
                     0.00976004 * pow(x, 4.) +
                     0.00107857 * pow(x, 5.);
    else if (x > 1)
        return (
            pow(x, 4.) +
            pow(x, 3.) * 8.5733287401 +
            pow(x, 2.) * 18.0590169730 +
            x * 8.6347608925 +
            0.2677737343) / (
            pow(x, 4.) + 
            pow(x, 3.) * 9.5733223454 +
            pow(x, 2.) * 25.6329561486 +
            x * 21.0996530827 +
            3.9584969228) / (
            x * exp(x));
    else
        return INFINITY;
}

double e2(double x) {
    if (x == 0.)
        return 1.;
    return exp(-x) - x * e1(x);
}

double hjerting(double u, double a) {
    int i;
    double hjertwing, hjertcore;

    u = fabs(u) + 1e-50;

    for (i = 0; i < (int) (sizeof(hjertu) / sizeof(double)) && u > hjertu[i]; i++);

    // printf("i: %d\tu: %.3lf\thjertu[i]: %.3lf\n",i, u, hjertu[i]);

    // printf("i: %i\n",i);
    if (i >= (int) (sizeof(hjertu) / sizeof(double))) {
        // printf("wing regime\n");
        hjertwing = (0.56419 * pow(u, -2.) + 0.846 * pow(u, -4.)) * a +
                (-0.56 * pow(u, -4.)) * pow(a, 3.);
        return hjertwing;
    }

    u = (u - hjertu[i-1]) / (hjertu[i] - hjertu[i-1]);

    hjertcore = ((hjert[0][i] - hjert[0][i-1]) * u + hjert[0][i-1]) +
                ((hjert[1][i] - hjert[1][i-1]) * u + hjert[1][i-1]) * a +
                ((hjert[2][i] - hjert[2][i-1]) * u + hjert[2][i-1]) * pow(a, 2.) +
                ((hjert[3][i] - hjert[3][i-1]) * u + hjert[3][i-1]) * pow(a, 3.) +
                ((hjert[4][i] - hjert[4][i-1]) * u + hjert[4][i-1]) * pow(a, 4.);

    return hjertcore;
}

static const void *dummy_ref[] = {dummy_ref, &c, &G, &h, &k, &mu, &me, &qe,
    &qeesu, &tsol, &Ihminus, &a0, &R, &fosc, &c4, &c6, &gammanat,
    speciesnames, A, I, amu, thetas, u0, u1, contjumps};