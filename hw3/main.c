#include <math.h>
#include <float.h>
#include <regex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "parse.h"
#include "phys.h"

/* EXTERNAL GLOBAL VARIABLES */
extern char *optarg;

/* GLOBAL VARIABLES */
char *f = "kurucz.txt";
double logg = 4.44;
double teff = 5777.;
double mtl[NSPECIES] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
char verbose = 0;
double pfnlines[NSPECIES][10][2][2];
double sahaphicurr[NSPECIES];
double sahaphihminus;

/* FUNCTION PROTOTYPES */
int parse(int argc, char **argv);
void init_pfn(void);
void init_mtl(void);
double pfn_interp(int species, int state, double temp);
void sahaphi(double temp);
double peguess(double pg);
double pe(double pg, double pei);
double pg(double mcol);
double *popdensities(double pg, double pe, double temp);
double *computemodel(double mcol, double temp);


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
    while((c = getopt(argc, argv, ":f:g:hmt:v")) != -1) {
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
                    -t\teffective temperature=5777.\n\
                    -v\tverbose=FALSE\n");
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

// make approx btwn each theta. at large t (small theta) extend linearly.
// at small t (large theta) use last value (u[j][9])
void init_pfn(void) {
    unsigned int j, t;
    double b0, m0, b1, m1;
    for(j = 0; j < NSPECIES; j++) {
        // printf("\nspecies %i\n",j);
        for(t = 0; t < sizeof(thetas)/sizeof(double) - 1; t++) {
            // printf("theta %i: %lf\n",t, thetas[t]);
            // printf("\nthetas[t+1]: %lf\tthetas[t]: %lf\n",thetas[t+1],thetas[t]);
            m0 = (u0[j][t+1] - u0[j][t]) / (thetas[t+1] - thetas[t]);
            b0 = u0[j][t] - m0 * thetas[t];

            // printf("u0[t+1]: %lf\tu0[t]: %lf\n",u0[j][t+1],u0[j][t]);
            // printf("m0: %lf\tb0: %lf\n",m0,b0);

            pfnlines[j][t][0][0] = m0;
            pfnlines[j][t][0][1] = b0;

            m1 = (u1[j][t+1] - u1[j][t]) / (thetas[t+1] - thetas[t]);
            b1 = u1[j][t] - m1 * thetas[t];

            // printf("u1[t+1]: %lf\tu1[t]: %lf\n",u1[j][t+1],u1[j][t]);
            // printf("m1: %lf\tb1: %lf\n",m1,b1);

            pfnlines[j][t][1][0] = m1;
            pfnlines[j][t][1][1] = b1;
        }

        pfnlines[j][9][0][0] = 0.;
        pfnlines[j][9][0][1] = u0[j][9];

        pfnlines[j][9][1][0] = 0.;
        pfnlines[j][9][1][1] = u1[j][9];

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
 * Yields the value of the partition function interpolated linearly on a log
 * scale. The returned value is the value of the partition function after
 * interpolation and exponentiation.
 */
double pfn_interp(int species, int state, double temp) {
    unsigned int t;
    double theta = 5040. / temp, m, b, u;
    
    // find which range the temperature belongs for determining pfn
    for(t = 0; t < sizeof(thetas)/sizeof(double); t++) {
        if(theta <= thetas[t])
            break;
    }
    
    // edge case where theta > 2.0 - use pfn for theta = 2.0
    if(t >= sizeof(thetas)/sizeof(double))
        return pow(10.,pfnlines[species][9][state][1]);

    // interpolate btwn the two nearest thetas and their pfn vals on log scale
    m = pfnlines[species][t][state][0];
    b = pfnlines[species][t][state][1];
    u = m * theta + b;
    return pow(10., u);
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
        sahaphicurr[j] = 0.6665 * pfn_interp(j, 1, temp)/pfn_interp(j, 0, temp) * pow(temp, 5./2.) *
            pow(10., -5040.*I[j]/temp);
    }
    sahaphihminus = 0.6665 * pfn_interp(j, 0, temp)/1.0 * pow(temp, 5./2.) *
            pow(10., -5040.*Ihminus/temp);
}

/**
 * Calculates a guess for the electron pressure Pe at this level given the
 * photosphere is only comprised of Hydrogen. This serves ONLY as an initial
 * value for subsequent iterations of the value of Pe.
 */
double peguess(double pg) {
    return -1. * sahaphicurr[0] + sqrt(pow(sahaphicurr[0],2.) + sahaphicurr[0] * pg);
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
 * Calculates the gas pressures Pg given the log-gravity and the mass column.
 * The calculation assumes a Kurucz atmosphere, where mass column is already
 * determined - no integration or iteration.
 */
double pg(double mcol) {
    return pow(10., logg) * mcol;
}

/**
 * Calculates the population densities of all neutral and ion species
 * considered in the model IN ADDITION TO H- ion and e-.
 */
double *popdensities(double pgf, double pef, double temp) {
    int j;
    double *pops, *ptr, asum, nh, nh1, nspeciestot, nspecies0, nspecies1;

    // Location where species populations are stored, including e- and H-
    pops = (double *) malloc((2 * NSPECIES + 2) * sizeof(double));
    ptr = pops;

    // electron number density
    *ptr = pef / (k * temp);    

    // H population statistics
    asum = 0.0;
    for (j = 0; j < NSPECIES; j++)
        asum += A[j];

    nh = (pgf - pef) / (k * temp) / asum;
    nh1 = nh / (pef/sahaphihminus + 1. + sahaphicurr[0]/pef);

    // H- number density
    ptr++;
    *ptr = nh1 * pef / sahaphihminus;

    // H I number density
    ptr++;
    *ptr = nh1;

    // H II number density
    ptr++;
    *ptr = nh1 * sahaphicurr[0] / pef;

    // population of all elements considered, excluding H-/HI/HII
    for (j = 1; j < NSPECIES; j++) {
        nspeciestot = nh * A[j];
        nspecies0 = nspeciestot / (1. + sahaphicurr[j]/pef);
        nspecies1 = nspeciestot - nspecies0;
        
        ptr++;
        *ptr = nspecies0;
        ptr++;
        *ptr = nspecies1;
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

    // Parse atmosphere file (kurucz model)
    nrows = parsekurucz(f, &databuf);

    // Initialize pfn and abundances via metallicity
    init_pfn();
    init_mtl();

    char *speciesnames[] = {
        "e-", "H-", "H I", "H II",
        "C I", "C II", "Na I", "Na II",
        "Mg I", "Mg II", "Si I", "Si II",
        "K I", "K II", "Ca I", "Ca II",
        "Cr I", "Cr II", "Fe I", "Fe II",
        "Ni I", "Ni II"     
    };

    // Calculate population densities for all levels in the Kurucz atmosphere
    for(i = 0; i < nrows; i++) {
        printf("\nLevel %i - %.5lf - %.2lf\n",i,*(databuf + 2 * i), *(databuf + 2 * i + 1));
        pops = computemodel(*(databuf + 2 * i), *(databuf + 2 * i + 1));
        printf("---- SPECIES POPULATIONS ----\nspecies\t\t#/cm^3\n");
        for(j = 0; j < 2 * NSPECIES + 2; j++) {
            printf("%s\t\t%.3le\n",speciesnames[j], *(pops + j));
        }
    }
}