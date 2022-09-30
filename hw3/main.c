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
double mtl[NSPECIES] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
long niters = 5;
char *f = "kurucz.txt";
char verbose = 0;
double sahaphicurr[NSPECIES];

int parse(int argc, char **argv);
double pg(double mcol);
void sahaphi(double temp);
double peguess(double pg);
double pe(double pg, double pei);
double *computemodel(double mcol, double temp);
void init_pfn(void);
void init_mtl(void);
double *popdensities(double pg, double pe, double temp);

/**
 * Calculates the gas pressures Pg given the log-gravity and the mass column.
 * The calculation assumes a Kurucz atmosphere, where mass column is already
 * determined - no integration or iteration.
 */
double pg(double mcol) {
    return pow(10., logg) * mcol;
}

/**
 * Calculates a guess for the electron pressure Pe at this level given the
 * photosphere is only comprised of Hydrogen. This serves ONLY as an initial
 * value for subsequent iterations of the value of Pe.
 */
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
    }
}

/**
 * Initializes abundance data based off of metallicity data.
 */
void init_mtl(void) {
    int i;
    for(i = 0; i < NSPECIES; i++) {
        if (mtl[i] != 0)
            A[i] = pow(10., mtl[i]) * A[i];
    }
}

/**
 * Calculates the Phi(T) expression in the Saha equation. The result is stored
 * in an array that holds all Phi(T) expressions for the current temperature.
 * The values in the array must be updated for every new level (temperature) in
 * the atmosphere for accurate calculations.
 */
void sahaphi(double temp) {
    int j;
    for(j = 0; j < NSPECIES; j++) {
        // uses theta=1.0 for all calcs currently
        sahaphicurr[j] = 0.6665 * u1[j][5]/u0[j][5] * pow(temp, 5./2.) *
            pow(10., -5040.*I[j]/temp);
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

    /* Iterates through all species to calculate their contribution to the
       numerator and denominator of the expression for Pe */
    for(j = 0; j < NSPECIES; j++) {
        factor = 1. / (pei / sahaphicurr[j] + 1.);
        numer += A[j] * factor;
        denom += A[j] * (1. + factor);
    }

    return pg * numer / denom;
}

/**
 * Calculates the population densities of all neutral and ion species
 * considered in the model IN ADDITION TO H- ion and e-.
 */
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

/**
 * Solves the atmosphere for a given level, determined by its mass column and
 * temperature. The gas pressure, electron pressure, and population densities
 * are calculated.
 */
double *computemodel(double mcol, double temp) {
    double pgi, pei, pef, precision = 1.e-7, *pops;
    int itercount = 0;

    // pre-calculate saha eqn Phi(T) for this level (temperature)
    sahaphi(temp);
    
    // calculate initial values for gas/electron pressure
    pgi = pg(mcol);
    pei = peguess(pgi);

    if (verbose)
        printf("pg=%le\nPe initl: %le\n",pgi, pei);

    // set initial guess to the working value of Pe
    pef = pei;
    
    // iteratively calculate Pe until it converges within desired precision
    do {
        pei = pef;
        pef = pe(pgi, pei);
        itercount++;
    } while(fabs(pei - pef)/pef > precision);

    if (verbose) {
        printf("Pe final: %le\n",pef);
        printf("after %i iterations, converged with at least precision %le\n",itercount,precision);
    }

    // determine population densities for each species desired
    pops = popdensities(pgi, pef, temp);

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
    while((c = getopt(argc, argv, ":f:g:hn:mt:v")) != -1) {
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
                    -m\tmetallicity={0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}\n\
                    -n\tnumber of iterations=1e3\n\
                    -t\teffective temperature=5777.\n");
                return 1;
            case 'm':
                // Parse metallicity (mtl) argument
                int mtlnum;
                double mtlval;

                do {
                    printf("Which metallicity would you like to update?\n\
                        0 - NONE\n\
                        1 - HYDROGEN\n\
                        2 - CARBON\n\
                        3 - SODIUM\n\
                        4 - MAGNESIUM\n\
                        5 - SILICON\n\
                        6 - POTASSIUM\n\
                        7 - CALCIUM\n\
                        8 - CHROMIUM\n\
                        9 - IRON\n\
                        10 - NICKEL\n");
                    scanf(" %i",&mtlnum);

                    if(mtlnum == 0)
                        break;
                    if(abs(mtlnum - 5) > 5) {
                        printf("Type an integer value from 0-10. You typed: %i\n",mtlnum);
                    } else {
                        printf("New metallicity: ");
                        scanf(" %le", &mtlval);
                        mtl[mtlnum - 1] = mtlval;
                    }
                } while (1);
                int i;
                printf("mtl: {%e",mtl[0]);
                for (i = 1; i < NSPECIES; i++) {
                    printf(", %e", mtl[i]);
                }
                printf("}\n");
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
            case 'v':
                verbose = 1;
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
    init_mtl();

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