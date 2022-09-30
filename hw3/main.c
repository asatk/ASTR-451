#include <math.h>
#include <float.h>
#include <regex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "parse.h"
#include "phys.h"

extern char *optarg;

double logg = 4.44;
double teff = 5777.;
double mtl = 0.;
long niters = 5;
char *f = "kurucz.txt";
double sahaphicurr[NSPECIES];

int parse(int argc, char **argv);
double pg(double mcol);
void sahaphi(double temp);
double peguess(double pg);
double pe(double pg, double pei);
double *computemodel(double mcol, double temp);
void init_pfn(void);
double *popdensities(double pg, double pe, double temp);

// MUST DO METALLICITY ADJUSTMENTS TO ABUNDANCES


// gas pressure determined from input parameters
double pg(double mcol) {
    return pow(10., logg) * mcol;
}

// guess of initial electron pressure assuming Hydrogen-only atmosphere
double peguess(double pg) {
    return -1. * sahaphicurr[0] + sqrt(pow(sahaphicurr[0],2.) + sahaphicurr[0] * pg);
}


// make approx btwn each theta. at large t (small theta) extend linearly.
// at small t (large theta) use last value (u[j][9])
void init_pfn(void) {
    int j;
    for(j = 0; j < NSPECIES; j++) {
        u0[j][5] = pow(10.,u0[j][5]);
        u1[j][5] = pow(10.,u1[j][5]);
        // printf("species %i: u0=%le\tu1=%le\n",j,u0[j],u1[j]);
    }
}

/**
 * 
 */
void sahaphi(double temp) {
    
    // double v0, v1, v2;
    int j;
    for(j = 0; j < NSPECIES; j++) {
        // v0 = 0.6665 * u1[j]/u0[j];
        // v1 = pow(temp, 5./2.);
        // v2 = pow(10., -5040.*I[j]/temp);
        // printf("v0=%le\tv1=%le\tv2=%le\tprod=%le\n",v0,v1,v2,v0*v1*v2); 
        // uses theta=1.0 for all calcs
        sahaphicurr[j] = 0.6665 * u1[j][5]/u0[j][5] * pow(temp, 5./2.) *
            pow(10., -5040.*I[j]/temp);
        // printf("species %i: sahaphi=%le\n",j,sahaphicurr[j]);
    }
}

/**
 * Calculates Pe (electron pressure) given a gas pressure and a guess at the
 * electron pressure. The final solved electron pressure is given in cgs units.
 * 
 * pre-req - saha phi run so that temp isnt important
 */
double pe(double pg, double pei) {
    int j;
    double numer = 0.0, denom = 0.0, factor;

    for(j = 0; j < NSPECIES; j++) {
        // printf("species %i: sahaphi=%le\n",j,sahaphicurr[j]);
        factor = 1. / (pei / sahaphicurr[j] + 1.);
        numer += A[j] * factor;
        denom += A[j] * (1. + factor);
    }

    return pg * numer / denom;
}

double *popdensities(double pg, double pe, double temp) {
    int j;
    double *pops, *ptr, asum, nh;

    // Location where species populations are stored, including e- and H-
    pops = (double *) malloc((2 * NSPECIES + 2) * sizeof(double));
    ptr = pops;

    // electron number density
    *ptr = pe / (k * temp);    

    // H- number density
    ptr++;
    *ptr = 1.0;

    // population of H
    asum = 0.0;
    for (j = 0; j < NSPECIES; j++)
        asum += A[j];

    nh = (pg - pe) / (k * temp) / asum;


    // population of all elements
    for (j = 0; j < NSPECIES; j++) {
        ptr++;
        *ptr = nh * A[j];
        ptr++;
        *ptr = nh * A[j];
    }

    

    return pops;


}

double *computemodel(double mcol, double temp) {
    double precision = 1.e-7;
    int itercount = 0;
    
    sahaphi(temp);
    
    double pgi = pg(mcol);

    printf("pg=%le\n",pgi);

    double pei = peguess(pgi);

    printf("Pe initl: %le\n",pei);
    
    // int i;
    /* Iteratively calculate Pg until its value converges to within 1% of the
       previous iteration's value */
    // for(i = 0; i < niters; i++) {
    //     // printf("electron pressure (Pe) guess: %le\n",pei);
    //     pei = pe(pgi, pei);
    // }
   
    double pef = pei;
    do {
        pei = pef;
        pef = pe(pgi, pei);
        itercount++;
    } while(fabs(pei - pef)/pef > precision);
    pei = pef;
    printf("Pe final: %le\n",pei);
    printf("after %i iterations, converged with at least precision %le\n",itercount,precision);

    double *pops;
    pops = popdensities(pgi, pei, temp);

    return pops;
}

/**
 * Parse the command-lne arguments passed in execution of the program.
 * 
 * argc     number of arguments, main()
 * argv     argument list from main()
 */
int parse(int argc, char **argv) {

    double val;
    int c;

    // Parse each character in the command-line argument string
    while((c = getopt(argc, argv, ":f:g:hn:m:t:")) != -1) {
        switch (c) {
            case 'f':
                // Parse filename (f) argument
                f = optarg;
                break;
            case 'g':
                // Parse gravity (logg) argument
                if ((val = parsenum(optarg)) != 0.0)
                    logg = val;
                printf("logg: %lf\n",logg);
                break;
            case 'h':
                // Print help menu
                printf("The following options and their default values are listed below:\n\
                    -f\tfilename=\"kurucz.txt\"\n\
                    -g\tlog gravity=4.44\n\
                    -h\thelp menu\n\
                    -m\tmetallicity=0.\n\
                    -n\tnumber of iterations=1e3\n\
                    -t\teffective temperature=5777.\n");
                return 1;
            case 'm':
                // Parse metallicity (mtl) argument
                if ((val = parsenum(optarg)) != 0)
                    mtl = val;
                printf("mtl: %lf\n",mtl);
                break;
            case 'n':
                // Parse number of iterations (niters) argument
                if ((val = parsenum(optarg)) != 0.0)
                    niters = (long) val;
                printf("niters: %ld\n",niters);
                break;
            case 't':
                // Parse temperature (teff) argument
                if ((val = parsenum(optarg)) != 0)
                    teff = val;
                printf("teff: %lf\n",teff);
                break;
            case ':':
                printf("option needs an argument\n");
                break;
            case '?':
                printf("unknown option\n");
                break;
        }
    }
    return 0;
}

/**
 * Entry point for the program. Computes spectrum for a star with parameters
 * specified in the command line for execution. Returns 0 upon successful exit.
 * 
 * argc     number of arguments passed in through command line
 * argv     list of arguments passed in through command line 
 */
int main(int argc, char **argv) {

    int err, nrows, i, j;
    double *databuf, *pops;

    // Parse the command-line arguments
    if ((err = parse(argc, argv)) != 0) {
        if (err != 1)
            printf("error parsing\n");
        return err;
    }

    nrows = parsekurucz(f, &databuf);
    (void)nrows;

    init_pfn();

    //replace nrows w 1
    for(i = 44; i < 45; i++) {
        printf("\nentry %i - %lf - %lf\n",i,*(databuf + 2 * i), *(databuf + 2 * i + 1));
        pops = computemodel(*(databuf + 2 * i), *(databuf + 2 * i + 1));
        for(j = 0; j < 2 * NSPECIES + 2; j++) {
            (void) pops;
            printf("pop of item %i: %le\n",j, *(pops + j));
        }
    }
}