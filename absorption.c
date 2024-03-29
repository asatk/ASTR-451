/**
 * 
 * Author: Anthony Atkinson
 * Date: 2022-13-10
*/

#include <math.h>
#include <stdio.h>

#include "phys.h"

/**
 * Calculates the excitation energy of the Hydrogen atom for a given principal
 * quantum number n. The result returned is in eV.
*/
static double chi(int n) {
    return I[HYDROGEN] * (1. - 1. / pow(n, 2.));
}

/**
 * Calculates the Gaunt Factor for Bound-Free Hydrogen absoprtion coefficient.
*/
static double g_bf(int n, double lambdaR) {
    // Eq. 8.5 (Gray 3ed)
    return 1 - 0.3456 * pow(lambdaR, -1./3.) * (lambdaR / pow(n, 2) - 0.5);
}

/**
 * Calculates the continuum absorption due to Hydrogen Bound-Free processes.
*/
static double k_hbf(double temp, double lambda) {
    int n, n0;
    double theta, sum;
    
    for (n0 = 0; n0 < (int) (sizeof(contjumps) / sizeof(double)); n0++) {
        if (lambda < contjumps[n0])
            break;
    }
    n0++;

    theta = 5040. / temp;

    sum = 0.;
    for (n = n0; n <= n0 + 2; n++) {
        sum += g_bf(n, lambda * R) * pow(10., -theta * chi(n)) / pow(n, 3.);
    }

    // Eq. 8.8 (Gray 3ed)
    return a0 * pow(lambda, 3.) * (sum + M_LOG10E / (2 * theta * I[HYDROGEN]) *
            (pow(10., -theta * chi(3)) - pow(10., -theta * I[HYDROGEN])));

}

/**
 * Calculates the continuum absorption due to Hydrogen Free-Free processes
*/
static double k_hff(double temp, double lambda) {
    double g_ff, theta;

    theta = 5040. / temp;

    // Eq. 8.6 (Gray 3ed)
    g_ff = 1 + 0.3456 * pow(lambda * R, -1./3.) * ((M_LOG10E * lambda) / (1.2398e4 * theta) + 0.5);

    // Eq. 8.10 (Gray 3ed)
    return a0 * pow(lambda, 3.) * g_ff * M_LOG10E / (2 * theta * I[HYDROGEN]) * pow(10., -theta * I[HYDROGEN]);
}

/**
 * H- Bound-Free absorption in cm^2 per neutral H atom.
 */
static double k_hminusbf(double temp, double pef, double lambda) {
    double a_bf, theta;

    theta = 5040. / temp;

    // Eq. 8.11 (Gray 3ed)
    a_bf =  1.99654 -                       // a0
            1.18267e-05 * lambda +          // a1 * lam
            2.64243e-06 * pow(lambda, 2.) - // a2 * lam^2
            4.40524e-10 * pow(lambda, 3.) + // a3 * lam^3
            3.23992e-14 * pow(lambda, 4.) - // a4 * lam^4
            1.39568e-18 * pow(lambda, 5.) + // a5 * lam^5
            2.78701e-23 * pow(lambda, 6.);  // a6 * lam^6

    // Eq 8.12 (Gray 3ed)
    return 4.158e-10 * 1.e-18 * a_bf * pef * pow(theta, 5./2.) * pow(10., Ihminus * theta);
}

/**
 * H- Free-Free absorption in cm^2 per neutral H atom.
 */
static double k_hminusff(double temp, double pef, double lambda) {
    double f0, f1, f2, loglambda, logtheta;

    loglambda = log10(lambda);
    logtheta = log10(5040. / temp);

    f0 =    -2.2763 - 
            1.6850 * loglambda + 
            0.76661 * pow(loglambda, 2.) -
            0.053346 * pow(loglambda, 3.);

    f1 =    15.2827 -
            9.2846 * loglambda +
            1.99381 * pow(loglambda, 2.) -
            0.142631 * pow(loglambda, 3.);

    f2 =    -197.789 +
            190.266 * loglambda -
            67.9775 * pow(loglambda, 2.) +
            10.6913 * pow(loglambda, 3.) -
            0.625151 * pow(loglambda, 4.);

    // Eq 8.13 (Gray 3ed)
    return 1.e-26 * pef *
            pow(10., f0 + f1 * logtheta + f2 * pow(logtheta, 2.));

}

/**
 * Calculates contiuum absorption due to e- Free-Free processes
*/
static double k_e(double pef, double pgf) {
    int j;
    double asum;

    asum = 0.0;
    for (j = 0; j < NSPECIES; j++)
        asum += A[j];

    // Eq 8.17 (Gray 3ed)
    return 0.6648e-24 * pef / (pgf - pef) * asum;
}

/**
 * Calculates the total continuum mass absorption coefficient in cm^2 g^-1 per H
*/
double k_total(double temp, double pef, double pgf, double lambda, double sahaphiH) {
    int j;
    double theta, sef, phpf, k_tot, msum;

    theta = 5040. / temp;

    sef = 1 - pow(10., -1.2398e4 * theta / lambda);
    phpf = 1. / (1. + sahaphiH / pef);

    // Eq. 8.18 (Gray 3ed)
    k_tot = ((k_hbf(temp, lambda) + k_hff(temp, lambda) +
        k_hminusbf(temp, pef, lambda)) * sef + k_hminusff(temp, pef, lambda)) *
        phpf + k_e(pef, pgf);

    msum = 0.0;
    for (j = 0; j < NSPECIES; j++)
        msum += A[j] * mu * amu[j];

    // Eq. 8.19 (Gray 3ed)
    return k_tot / msum;
}

/**
 * Calculates the total line mass absorption coefficient in cm^2 g^-1 per H
*/
double ell_total(double temp, double pef, double pgf, double lambda,
    double linecenter, int species, double sahaboltzelement) {
    
    int j;
    double a, dlam, dlamdopp, gamma4, gamma6, gamma_tot, msum, sef, u;

    // Sum of masses of each element in the star relative to their abundances
    msum = 0.0;
    for (j = 0; j < NSPECIES; j++)
        msum += A[j] * mu * amu[j];

    // Stimulated-emission factor
    sef = (1 - pow(10., -1.2398e4 * 5040. / temp / lambda));

    // Collision broadening - quadratic Stark
    gamma4 = pow(10., 19. + 2./3. * c4 + log10(pef) - 5./6. * log10(temp));
    // Collision broadening - van der Waals
    gamma6 = pow(10., 20. + 0.4 * c6 + log10(pgf) - 0.7 * log10(temp));
    // Total Lorentz width/damping
    gamma_tot = gamma4 + gamma6 + pow(10., gammanat);
    // Doppler width (Gaussian), Eq. 11.39
    dlamdopp = 4.301e-7 * (linecenter) * pow(temp / (amu[species]), 1. / 2.);
    // Distance from line center
    dlam = lambda - linecenter;

    // Hjerting function parameters
    // Eq. 11.47/11.55 (Gray 3ed)
    a = 2.65e-20 * gamma_tot * pow(linecenter, 2.) / dlamdopp;
    u = dlam / dlamdopp;

    // Eq. 11.54 (Gray 3ed)   
    return 4.995e-21 * (hjerting(u, a) * pow(linecenter, 2.) * A[species] *
        fosc * sahaboltzelement) / (dlamdopp * msum) * sef;
}


// Prevent "unused function" and "unused variable" warnings.
static const void *dummy_ref[] = {amu, &fosc, &tsol, &qe, &qeesu, &me, &mu, &G, &h, &c,
    &k, &gammanat, &c4, &c6, u0, u1, thetas, speciesnames, hjertu, hjert, dummy_ref};