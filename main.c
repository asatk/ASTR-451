/**
 * Computes a model photosphere of a star with specified parameters. Model
 * is based off of a Kurucz atmosphere, with gas pressure computed easily from
 * the mass column. The populations present in the atmosphere are determined
 * by the temperature and mass column.
 * 
 * Author: Anthony Atkinson
 * Date: 2021-09-30
 */

#include <math.h>
#include <float.h>
#include <regex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// #include "absorption.h"
#include "parse.h"
#include "phys.h"
#include "absorption.h"

/* EXTERNAL GLOBAL VARIABLES */
extern char *optarg;

/* GLOBAL VARIABLES */
char *f = "kurucz.txt";
char *outfilenamepops = "pops.dat";
char *outfilenameflux = "flux_lowg.dat";
double logg = 4.44;
double teff = 5777.;
double mtl[NSPECIES] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
char logscale = 0;
int plotmode = 0;
char verbose = 0;
double pfnlines[NSPECIES][10][2][2]; // [species][theta index][state][m=0,b=1]
double sahaphicurr[NSPECIES];
double sahaphihminus;
double linemin = 5885.;
double linemax = 5895.;
double linecnt = 5890.;
double lineinc = 0.01;
double *tauccurr;
double *taulcurr;

/* FUNCTION PROTOTYPES */
int parse(int argc, char **argv);
void init_hfn(void);
void init_pfn(void);
void init_mtl(void);
double pfn_interp(int species, int state, double temp);
void sahaphi(double temp);
double peguess(double pg);
double pe(double pg, double pei);
double pg(double mcol);
double *popdensities(double pg, double pe, double temp);
double *computemodel(double mcol, double dmcol, double temp, double *flux);
void printlevel(double *pops, int level, FILE *outfile);
double *lineflux(double *flux, double mcol, double temp, double pef, double pgf);
void plotpops();


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
    while((c = getopt(argc, argv, ":f:g:hmopt:v")) != -1) {
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
                    -m\tmetallicity={0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}\n\
                    -o\tlogscale=FALSE\n\
                    -p\tplotmode=0\n\
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
                        2 - HELIUM\n\
                        3 - CARBON\n\
                        4 - SODIUM\n\
                        5 - MAGNESIUM\n\
                        6 - SILICON\n\
                        7 - POTASSIUM\n\
                        8 - CALCIUM\n\
                        9 - CHROMIUM\n\
                        10 - IRON\n\
                        11 - NICKEL\n");
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
            case 'o':
                logscale = 1;
                break;
            case 'p':
                printf("Which plotting mode?\n\
                    0 - OFF\n\
                    1 - PNG\n\
                    2 - PDF\n");
                scanf(" %i",&plotmode);
                if (plotmode == 2)
                    printf("Plotting PDF only works if the program 'ps2pdf' is installed\n");
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
 * Approximate the partition function by linearly interpolating between y=log u
 * and x=theta values. At large theta, want asymptotic behavior thus for any
 * theta greater than 2.0, its pfn value is that at 2.0. At small theta, want
 * increasing behavior as temperature grows, thus its pfn value is scaled
 * linearly with previous interpolation. Stores interpolated parameters in the
 * 4D array pfnlines where indices are [species][theta index][state][m=0,b=1].
 */
void init_pfn(void) {
    unsigned int j, t;
    double b0, m0, b1, m1;
    // calculate interpolation parameters for each species
    for(j = 0; j < NSPECIES; j++) {
        /* determine slope and intercept of line connecting consecutive points
           from the tabulated values in Table D.2 (Gray 3ed) */
        for(t = 0; t < sizeof(thetas)/sizeof(double) - 1; t++) {
            // neutral species pfn
            m0 = (u0[j][t+1] - u0[j][t]) / (thetas[t+1] - thetas[t]);
            b0 = u0[j][t] - m0 * thetas[t];

            pfnlines[j][t][0][0] = m0;
            pfnlines[j][t][0][1] = b0;

            // singly-ionized species pfn
            m1 = (u1[j][t+1] - u1[j][t]) / (thetas[t+1] - thetas[t]);
            b1 = u1[j][t] - m1 * thetas[t];

            pfnlines[j][t][1][0] = m1;
            pfnlines[j][t][1][1] = b1;
        }

        /* Boundary case: theta > 2.0. want asymptotic behavior for pfn;
           otherwise, for sufficiently low T (high theta), pfn is negative! */
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
 * Yields the value of the partition function interpolated linearly on a log u
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
    sahaphihminus = 0.6665 * pfn_interp(HYDROGEN, 0, temp)/1.0 * pow(temp, 5./2.) *
            pow(10., -5040.*Ihminus/temp);
}

/**
 * Calculates a guess for the electron pressure Pe at this level given the
 * photosphere is only comprised of Hydrogen. This serves ONLY as an initial
 * value for subsequent iterations of the value of Pe. Returns the guess.
 */
double peguess(double pg) {
    return -1. * sahaphicurr[HYDROGEN] + sqrt(pow(sahaphicurr[HYDROGEN],2.) + sahaphicurr[HYDROGEN] * pg);
}

/**
 * Calculates Pe (electron pressure) given a gas pressure and a guess at the
 * electron pressure. The final returned electron pressure is given in cgs
 * units.
 * Pre-condition: sahaphi was run at the given temperature level.
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
 * determined - no integration or iteration. Returns the calculated pg.
 */
double pg(double mcol) {
    return pow(10., logg) * mcol;
}

/**
 * Calculates the population densities of all neutral and ion species
 * considered in the model IN ADDITION TO H- ion and e-.
 * 
 * Returns the number densities for all desired species as a pointer to a
 * memory location with 2 * NSPECIES + 2 double elements. Each element contains
 * a corresponding number density: e-, H-, H I, H II, C I, C II, ...Ni I, Ni II
 */
double *popdensities(double pgf, double pef, double temp) {
    int j;
    double *pops, *ptr, asum, nh, nh1, nspeciestot, nspecies0, nspecies1;

    // Location where species populations are stored, including e- and H-
    pops = (double *) malloc((2 * NSPECIES + 2) * sizeof(double));
    ptr = pops;

    // electron number densityz
    *ptr = pef / (k * temp);    

    // H population statistics
    asum = 0.0;
    for (j = 0; j < NSPECIES; j++)
        asum += A[j];

    nh = (pgf - pef) / (k * temp) / asum;
    nh1 = nh / (pef/sahaphihminus + 1. + sahaphicurr[HYDROGEN]/pef);

    // H- number density
    ptr++;
    *ptr = nh1 * pef / sahaphihminus;

    // H I number density
    ptr++;
    *ptr = nh1;

    // H II number density
    ptr++;
    *ptr = nh1 * sahaphicurr[HYDROGEN] / pef;

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
 * 
 * Returns the number densities for all desired species as a pointer to a
 * memory location with 2 * NSPECIES + 2 double elements. Each element contains
 * a corresponding number density: e-, H-, H I, H II, C I, C II, ...Ni I, Ni II
 */
double *computemodel(double mcol, double dmcol, double temp, double *flux) {
    double pgf, pei, pef, precision = 1.e-10, *pops;
    int itercount = 0;

    // pre-calculate saha eqn Phi(T) for this level (temperature)
    sahaphi(temp);
    
    // calculate initial values for gas/electron pressure
    pgf = pg(mcol);
    pei = peguess(pgf);

    if (verbose)
        printf("pg=%le\nPe initl: %le\n",pgf, pei);

    // set initial guess to the working value of Pe
    pef = pei;
    
    // iteratively calculate Pe until it converges within desired precision
    do {
        pei = pef;
        pef = pe(pgf, pei);
        itercount++;
    } while(fabs(pei - pef)/pef > precision);

    if (verbose) {
        printf("Pe final: %le\n",pef);
        printf("after %i iterations, converged with at least precision %le\n",itercount,precision);
    }

    // determine population densities for each species desired
    pops = popdensities(pgf, pef, temp);

    lineflux(flux, dmcol, temp, pef, pgf);

    if (verbose)
        printf("continuum absorption coefficient for 5890A: %.3le\n", k_total(temp, pef, pgf, 5890., sahaphicurr[HYDROGEN]));

    return pops;
}

void printlevel(double *pops, int level, FILE *outfile) {
    int j;
    char *temp_name;

    if (verbose)
        printf("---- SPECIES POPULATIONS ----\n%-*s#/cm^3\n",10,"species");
    
    // H- to e-
    if (verbose)
        printf("%-*s%.3le\n",10,speciesnames[0], *(pops + 1) / *(pops + 0));
    temp_name = (char *) malloc(strlen(speciesnames[0]) + 2);
    temp_name = strcpy(temp_name, speciesnames[0]);
    temp_name = strcat(temp_name,"\'");
    fprintf(outfile, "%-*i\'%-*s%-*le\n",8,level,10,temp_name,12, *(pops + 1) / *(pops + 0));

    // print each species and its density to outfile and stdin if verbose
    for(j = 0; j < 2 * NSPECIES + 2; j++) {
        if (verbose)
            printf("%-*s%.3le\n",10,speciesnames[j + 1], *(pops + j));
        temp_name = (char *) malloc(strlen(speciesnames[j + 1]) + 2);
        temp_name = strcpy(temp_name, speciesnames[j + 1]);
        temp_name = strcat(temp_name,"\'");
        fprintf(outfile, "%-*i\'%-*s%-*le\n",8,level,10,temp_name,12,*(pops + j));
    }
}

double *lineflux(double *flux, double dmcol, double temp, double pef, double pgf) {
    int i, nsteps;
    double sahaboltzelement,ell, kappa, lambda, *ptrflux, *ptrtauc, *ptrtaul, tauc, taul,
        dtauc, dtaul, fluxc, fluxl;

    printf("TEMP: %-*.3lePE: %-*.3lePG: %-*.3le\n",12,temp,12,pef,12,pgf);

    nsteps = (linemax - linemin) / lineinc;

    /* Calculate the combined Saha-Boltzmann equation for the number density of
        atoms in the desired excitation and ionization state */
    sahaboltzelement = 2. * pow(10., -5040 / temp * 0.) / 
        pfn_interp(SODIUM, 0, temp) / (1. + sahaphicurr[SODIUM] / pef);

    ptrtauc = tauccurr;
    ptrtaul = taulcurr;
    ptrflux = flux;

    // Iterate through every wavelength in the desired range around the line
    for (i = 0; i <= nsteps; i++) {

        lambda = linemin + i * lineinc;
        
        // Calculate mass absorption coefficients
        kappa = k_total(temp, pef, pgf, lambda, sahaphicurr[HYDROGEN]);
        ell = ell_total(temp, pef, pgf, lambda, linecnt, SODIUM, sahaboltzelement);

        // Calculate the increment in optical depth for this level
        dtauc = kappa * dmcol;
        dtaul = (kappa + ell) * dmcol;

        // Increment optical depth from previous level
        *ptrtauc += dtauc;
        *ptrtaul += dtaul;

        tauc = *ptrtauc;
        taul = *ptrtaul;

        // Calculate fluxes at this level for a lambda - Eq. 7.15 (Gray 3ed)
        fluxc = 2 * M_PI * planck(temp, lambda * 1e-8) * e2(tauc) * dtauc;
        fluxl = 2 * M_PI * planck(temp, lambda * 1e-8) * e2(taul) * dtaul;
        
        *ptrflux += fluxl;
        ptrflux++;
        *ptrflux += fluxc;

        ptrtauc++;
        ptrtaul++;
        ptrflux++;
    }

    return flux;
}

/**
 * Plot the number densities from the output file pops.dat using gnuplot.
 */
void plotpops() {
    int j;
    char str[300], fmt_str[300], title[10], titlej[10];

    // modifications for logscale in y
    if (logscale) {
        strcpy(fmt_str, "gnuplot -e \"set logscale y 10; ");
        strcpy(title, "log ");
    } else {
        strcpy(fmt_str, "gnuplot -e \"");
        strcpy(title, "");
    }
    
    // gnuplot command template
    if (plotmode == 1)
        strcat(fmt_str,"set terminal png size 500,500; set output '%s.png'; set xlabel 'Photosphere Level'; set ylabel 'Number Density (%s per cm^{-3})'; plot '<(cat pops.dat | grep \\\"%s''\\\" | cat)' using 4 title '%s'\"");
    else if (plotmode == 2)
        strcat(fmt_str,"set terminal postscript; set output '| ps2pdf - \\\"%s.pdf\\\"'; set xlabel 'Photosphere Level'; set ylabel 'Number Density (%s per cm^{-3})'; plot '<(cat pops.dat | grep \\\"%s''\\\" | cat)' using 4 title '%s'\"");
    
    // H- to e-
    strcpy(titlej, title);
    strcat(titlej, "H-e- ");
    snprintf(str, 300, fmt_str, titlej, speciesnames[0], speciesnames[0], titlej);
    system(str);

    // plot each species
    for(j = 1; j < 2 * NSPECIES + 3; j++) {
        strcpy(titlej, title);
        strcat(titlej, speciesnames[j]);
        
        // complete template command
        snprintf(str, 300, fmt_str, titlej, speciesnames[j], speciesnames[j], titlej);
        
        // run command
        system(str);
    }
}

/**
 * Entry point for the program. Computes spectrum for a star with parameters
 * specified in the command line for execution. Returns 0 upon successful exit.
 * 
 * argc     number of arguments passed in through command line
 * argv     list of arguments passed in through command line 
 */
int main(int argc, char **argv) {

    int err, nrows, nsteps, i, j;
    double *databuf, *pops, *flux, *ptr, *ptrl, *ptrc, dmcol;
    FILE *outfilepops, *outfileflux;

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

    // Open output file to write
    outfilepops = fopen(outfilenamepops,"w");
    fprintf(outfilepops, "#%-7s%-8s%-8s\n#%-7s%-10s%-8s\n","level","species","density","(#)","","(cm^-3)");

    nsteps = (linemax - linemin) / lineinc;
    
    // Calculate population densities  and flux for all levels in Kurucz atmos
    flux = (double *) malloc(2 * (nsteps + 1) * sizeof(double));
    ptr = flux;
    for(j = 0; j < 2 * (nsteps + 1); j++) {
        *ptr = 0.;
        ptr++;
    }

    // Initialize optical depth arrays for every wavelength to 0
    tauccurr = (double *) malloc((nsteps + 1) * sizeof(double));
    taulcurr = (double *) malloc((nsteps + 1) * sizeof(double));
    
    ptrl = taulcurr;
    ptrc = tauccurr;
    for(j = 0; j < nsteps + 1; j++) {
        *ptrl = 0.;
        *ptrc = 0.;
        ptrl++;
        ptrc++;
    }

    /* Iterate through each depth and compute the number density of all relevant
        atoms and the flux spectrum for the desired line */
    for(i = 0; i < nrows; i++) {
        dmcol = (i == 0) ? *(databuf) : *(databuf + 2 * i) - *(databuf + 2 * (i - 1));
        if (verbose)
            printf("\nLevel %i - %.5le - %.2le\n",i,*(databuf + 2 * i), *(databuf + 2 * i + 1));
        pops = computemodel(*(databuf + 2 * i), dmcol, *(databuf + 2 * i + 1), flux);
        
        printlevel(pops, i, outfilepops);
    }

    fclose(outfilepops);

    if (plotmode != 0)
        plotpops();

    outfileflux = fopen(outfilenameflux, "w");

    // Output the total line and continuum fluxes at every level
    for(j = 0; j < 2 * (nsteps + 1); j+=2) {
        // printf("fluxl: %-*.3lefluxc: %-*.3lel/c %-*.3le\n",12,*(flux + j), 12, *(flux + j + 1), 12, *(flux + j) / *(flux + j + 1));
        fprintf(outfileflux, "%-*.4le%-*.4le\n",15,linemin + j * lineinc,15,*(flux + j) / *(flux + j + 1));
    }

    fclose(outfileflux);
}

// Prevent "unused function" and "unused variable" warnings.
static const void *dummy_ref[] = {amu, contjumps, &fosc, &tsol, &qe, &qeesu,
    &me, &mu, &G, &h, &c, &a0, &R, &gammanat, &c4, &c6, hjertu, hjert,
    dummy_ref};