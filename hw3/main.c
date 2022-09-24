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
double mtl = 0;
long nsteps = 1e3;
char *f = "kurucz.txt";

int parse(int argc, char **argv);
double tstar(double tau0);
double pg(double pgi, double pei, double tau0);
double sahaphi(double temp, int species);
double pe(double pg, double pei, double temp);
double kappa(double pe, double temp);
double computemodel(double mcol, double temp);

//sum of all absorption methods
double kappa(double pe, double temp) {
    double kappa0;
    
    return kappa0;
}

// temperature of star determined by grey atmosphere/solar model
double tstar(double tau0) {
    return teff * pow(3./4. * (tau0 + 2./3.), 1./4.);
}

// return optical depth from temperature based on grey atm model
double tau(double temp) {
    return pow(temp / teff, 4.) * 4./3. - 2./3.;
}

// gas pressure determined iteratively assuming hydrostatic equilibrium
double pg(double pgi, double pei, double tau0) {
    int i;
    double pgf, pintg, val, t0, dlogt0, pe0, temp0, kappa0;

    temp0 = tstar(tau0);
    pe0 = pe(pgi, pei, temp0);
    kappa0 = kappa(pe0, temp0);

    // double precision = 1e-7; //fraction precision to which the integral must converge
    pintg = 0.0;
    dlogt0 = (DBL_MAX - tau0)/nsteps;

    for(i = 0; i < nsteps; i++) {
        val = (pow(t0, 1./2.) * pow(pgi, 1./2.) / (kappa0 * M_LOG10E)) * dlogt0;
        t0 += pow(10., dlogt0);
        pintg += val;
    }

    printf("pintg: %lf\n", pintg);

    pgf = pow(3./2. * pow(10., logg) * pintg, 2./3.);

    return pgf;
}

/**
 * Returns the electron pressure-independent factor of the Saha equation
 * corresponding to a given temperature and ion species. Uses cgs units.
 */
double sahaphi(double temp, int species) {
    return 0.6665 * u1[species]/u0[species] * pow(temp, 5./2) * pow(10., -5040.*I[species]/temp);
}

/**
 * Calculates Pe (electron pressure) given a gas pressure, an initial guess at
 * electron pressure, and temperature. The final solved electron pressure is
 * given in cgs units.
 */
double pe(double pg, double pei, double temp) {
    unsigned int j;
    double pef, numer = 0.0, denom = 0.0, factor;

    for(j = 0; j < sizeof(enum species)/sizeof(int); j++) {
        factor = 1. / (pei/sahaphi(temp, j) + 1.);
        numer += A[j] * factor;
        denom += A[j] * (1. + factor);
    }
    pef = numer/denom;
    return pef;
}

double computemodel(double mcol, double temp) {
    double precision = 0.01;
    double pgi = mcol * pow(10., logg);
    double tau0 = tau(temp);
    double pgf;
    int i = 0;
    do {
        printf("iteration %i\n\tpgi: %lf\n",i,pgi);
        pgf = pg(pgi, tau0);
        printf("pgf: %lf",pgf);
    } while((pgi - pgf)/pgf > precision);
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
    while((c = getopt(argc, argv, ":g:n:m:t:")) != -1) {
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
            case 'n':
                // Parse number of iterations (nsteps) argument
                if ((val = parsenum(optarg)) != 0.0)
                    nsteps = (long) val;
                printf("nsteps: %ld\n",nsteps);
                break;
            case 'm':
                // Parse metallicity (mtl) argument
                if ((val = parsenum(optarg)) != 0)
                    mtl = val;
                printf("mtl: %lf\n",mtl);
                break;
            case 't':
                // Parse temperature (teff) argument
                if ((val = parsenum(optarg)) != 0)
                    teff = val;
                printf("teff: %lf\n",teff);
                break;
            default:
                printf("hello from default\n");
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

    int err, nrows, i;
    double *databuf;

    // Parse the command-line arguments
    if ((err = parse(argc, argv)) != 0) {
        printf("error parsing\n");
        return err;
    }

    nrows = parsekurucz(f, &databuf);

    //replace nrwos w 1
    for(i = 0; i < 1; i++) {
        printf("entry %i - %lf - %lf\n",i,*(databuf + 2 * i), *(databuf + 2 * i + 1));
        computemodel(*(databuf + 2 * i), *(databuf + 2 * i + 1));
    }

    

}