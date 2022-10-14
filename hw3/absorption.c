/**
 * 
 * Author: Anthony Atkinson
 * Date: 2022-13-10
*/

#include <math.h>

#include "phys.h"

a0 = 1.0449e-26;        // for lambda in angstroms
R = 1.0968e-3;          // for lambda in angstroms
loge = 0.43429;

/**
 * 
 * Continuum absorption: kappa
 * srcs:
 *  - H I bf (SEF) (PHPF)
 *  - H I ff (SEF) (PHPF)
 *  - H- bf (SEF) (PHPF)
 *  - H2+ (SEF) (PHPF)
 *  - H- ff (PHPF)
 *  - He- ff
 *  - e-
 *  - metals
 * 
 * SEF = needs stim emission correction: 1/(1 + sahaphi/pe)
 * HPF = needs per H ptcl correction 1-10^(-X(lambda)*theta)
 * 
 */

double chi(int n) {
    return I[HYDROGEN] * (1. - 1. / pow(n, 2.));
}

double g_bf(int n, double lambdaR) {
    // double R = 1.0968e-3;          // for lambda in angstroms
    // lambdaR is product

    return 1 - 0.3456 * pow(lambdaR, -1./3.) * (lambdaR / pow(n, 2) - 0.5);
}

double k_hbf(double temp, double lambda) {

    int n, n0;
    double a0, R, theta, sum;

    theta = 5040. / temp;

    // a0 = 1.0449e-26;        // for lambda in angstroms
    // R = 1.0968e-3;          // for lambda in angstroms

    //log e = 0.43429

    sum = 0.;
    for (n = n0; n < n0 + 3; n++) {
        sum += g_bf(n, lambda * R) * pow(10., theta * chi(n)) / pow(n, 3.);
    }

    return a0 * pow(lambda, 3.) * (sum + loge / (2 * theta * I[HYDROGEN]) * (pow(10., -theta * chi(3)) - pow(10., -theta * I[HYDROGEN])));

}

double k_hff(double temp, double lambda) {

    double a0, g_ff, theta;

    theta = 5040. / temp;

    g_ff = 1 + 0.3456 * pow(lambda * R, -1./3.) * ((loge * lambda) / (1.2398e4 * theta) + 0.5);

    return a0 * pow(lambda, 3.) * g_ff * loge / (2 * theta * I[HYDROGEN]) * pow(10., -theta * I[HYDROGEN]);
}


/**
 * H- Bound-Free absorption in cm^2 per neutral H atom.
 */
double k_hminusbf(double pef, double temp, double lambda) {
    double a_bf, theta;

    theta = 5040. / temp;

    a_bf =  1.99654 -                       // a0
            1.18267e-05 * lambda +          // a1 * lam
            2.64243e-06 * pow(lambda, 2.) - // a2 * lam^2
            4.40524e-10 * pow(lambda, 3.) + // a3 * lam^3
            3.23992e-14 * pow(lambda, 4.) - // a4 * lam^4
            1.39568e-18 * pow(lambda, 5.) + // a5 * lam^5
            2.78701e-23 * pow(lambda, 6.);  // a6 * lam^6

    return 4.158e-10 * a_bf * pef * pow(theta, 5./2.) * pow(10., Ihminus * theta);
}

/**
 * H- Free-Free absorption in cm^2 per neutral H atom.
 */
double k_hminusff(double pef, double temp, double lambda) {
    double f0, f1, f2, loglambda, logtheta;

    loglambda = log10(loglambda);
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

    return 1.e-26 * pef *
            pow(10., f0 + f1 * logtheta + f2 * pow(logtheta, 2.));

}

// double k_h2plus()

// double k_heminusff()

double k_e(double pef, double pgf) {
    int j;
    double asum;

    asum = 0.0;
    for (j = 0; j < NSPECIES; j++)
        asum += A[j];

    // this could also be just ratio of Ne- to NH depending on how I do this p 162 Gray
    return 0.6648e-24 * pef / (pgf - pef) * asum;
}

double k_metals() {
    return 0.;
}

double k_total(double pef, double pgf, double temp, double lambda, double sahaphiH) {
    int j;
    double theta, sef, phpf, k_tot, msum;


    // make this the param for abs
    theta = 5040. / temp;

    sef = 1 - pow(10., -1.2398e4 * theta / lambda);
    phpf = 1. / (1. + sahaphiH / pef);

    k_tot = ((k_hbf(temp, lambda) + k_hff(temp, lambda) + k_hminusbf(pef, temp, lambda)) * sef + k_hminusff(pef, temp, lambda)) * phpf + k_e(pef, pgf) + k_metals();

    msum = 0.0;
    for (j = 0; j < NSPECIES; j++)
        msum += A[j] * mu * amu[j];

    return k_tot / msum;
}

/**
 * Line absorption: ell
 * damping srcs:
 *  - natural broadening (L) gamma
 *  - pressure broadening (L) gamma4, gamma6
 *  - thermal/collisional braodening (G)
*/