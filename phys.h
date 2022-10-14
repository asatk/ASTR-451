/**
 * Definition of general physical constants as well as properties of a star
 * assumed true ubiquitously.
 * 
 * Author: Anthony Atkinson
 * Date: 2021-09-30
 */

#define NSPECIES 11

enum species {HYDROGEN, HELIUM, CARBON, SODIUM, MAGNESIUM, SILICON, POTASSIUM, CALCIUM,
    CHROMIUM, IRON, NICKEL};

static double c = 2.998e10;    // speed of light, cm s-1
static double G = 6.674e-8;    // gravitational constant, ergs cm g-2
static double k = 1.381e-16;   // boltzmann constant, J K-1
static double mu = 1.661e-24;  // atomic mass unit, g
static double me = 9.109e-28;  // electron mass, g
static double qe = 1.602e-19;  // electron charge, C 4.803e-10 e.s.u. in mks units p 233 Gray
static double tsol = 5.777e3;  // effective temperature of the sun (Sol)
static double Ihminus = 0.754; // ionization energy of H- ion
static double a0 = 1.0449e-26; // for lambda in angstroms
static double R = 1.0968e-3;   // for lambda in angstroms
static double loge = 0.43429;  // log10(e)
static double fosc = 0.6546;   // oscillator strength (f) used in natural broadening
static double c4 = -15.17;     // log C4 constant for 
static double c6 = -31.7;      // log C6 constant for
static double gammanat = 7.9;  // log gamma (damping) for natural broadening

/**
 * ABUNDANCES (# atom/# hydrogen, in the Sun)
 * per Table 16.3 (Gray, 3ed)
 */
static double A[] = {
    1.00e+0,    //H
    8.51e-2,    //He
    3.31e-4,    //C
    2.14e-6,    //Na
    3.80e-5,    //Mg
    3.35e-5,    //Si
    1.32e-7,    //K
    2.29e-6,    //Ca
    4.68e-7,    //Cr
    2.75e-5,    //Fe
    1.78e-6     //Ni
};

/**
 * IONIZATIONS (eV for 1st ionization)
 * per Table D.1 (Gray, 3ed)
 */
static double I[] = {
    13.598,     //H
    24.587,     //He
    11.260,     //C
    5.139,      //Na
    7.646,      //Mg
    8.152,      //Si
    4.341,      //K
    6.113,      //Ca
    6.767,      //Cr
    7.902,      //Fe
    7.640       //Ni
};

/**
 * ATOMIC WEIGHTS (amu)
 * per Table D.1 (Gray, 3ed)
 */
static double amu[] = {
    1.008,      //H
    4.003,      //He
    12.011,     //C
    22.990,     //Na
    24.305,     //Mg
    28.086,     //Si
    39.10,      //K
    40.08,      //Ca
    51.996,     //Cr
    55.847,     //Fe
    58.71       //Ni
};

/**
 * Theta values corresponding to the tabulated partition function argument
 */
static double thetas[10] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};

/**
 * PARTITION FUNCTIONS for NEUTRAL SPECIES
 * per Table D.2 (Gray, 3ed)
 */
static double u0[][10] = {
//   0.2    0.4    0.6    0.8    1.0    1.2    1.4    1.6    1.8    2.0
    {0.368, 0.303, 0.301, 0.301, 0.301, 0.301, 0.301, 0.301, 0.301, 0.301}, //H
    {0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000}, //He
    {1.163, 1.037, 0.994, 0.975, 0.964, 0.958, 0.954, 0.951, 0.950, 0.948}, //C
    {4.316, 1.043, 0.493, 0.357, 0.320, 0.309, 0.307, 0.306, 0.306, 0.306}, //Na
    {2.839, 0.478, 0.110, 0.027, 0.007, 0.002, 0.001, 0.001, 0.001, 0.000}, //Mg
    {1.521, 1.111, 1.030, 0.996, 0.976, 0.961, 0.949, 0.940, 0.932, 0.925}, //Si
    {4.647, 1.329, 0.642, 0.429, 0.351, 0.320, 0.308, 0.303, 0.302, 0.302}, //K
    {5.238, 1.332, 0.465, 0.181, 0.073, 0.028, 0.010, 0.003, 0.001, 0.000}, //Ca
    {4.284, 1.977, 1.380, 1.141, 1.022, 0.956, 0.917, 0.892, 0.875, 0.865}, //Cr
    {3.760, 2.049, 1.664, 1.519, 1.446, 1.402, 1.372, 1.350, 1.332, 1.317}, //Fe
    {2.779, 1.753, 1.577, 1.521, 1.490, 1.467, 1.447, 1.428, 1.410, 1.394}  //Ni
};

/**
 * PARTITION FUNCTIONS for SINGLY-IONIZIED SPECIES at theta=1.0
 * per Table D.2 (Gray, 3ed)
 */
static double u1[][10] = {
//   0.2    0.4    0.6    0.8    1.0    1.2    1.4    1.6    1.8    2.0
    {0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000}, //H+
    {0.301, 0.301, 0.301, 0.301, 0.301, 0.301, 0.301, 0.301, 0.301, 0.301}, //He+
    {0.853, 0.782, 0.775, 0.774, 0.773, 0.772, 0.771, 0.770, 0.769, 0.767}, //C+
    {0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000}, //Na+
    {0.537, 0.326, 0.304, 0.301, 0.301, 0.301, 0.301, 0.301, 0.301, 0.301}, //Mg+
    {0.900, 0.778, 0.764, 0.759, 0.755, 0.750, 0.746, 0.741, 0.736, 0.731}, //Si+
    {0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000}, //K+
    {0.825, 0.658, 0.483, 0.391, 0.344, 0.320, 0.309, 0.304, 0.302, 0.301}, //Ca+
    {1.981, 1.489, 1.125, 0.944, 0.856, 0.813, 0.793, 0.784, 0.781, 0.780}, //Cr+
    {2.307, 1.881, 1.749, 1.682, 1.638, 1.604, 1.575, 1.549, 1.525, 1.504}, //Fe+
    {1.659, 1.386, 1.215, 1.108, 1.037, 0.988, 0.953, 0.927, 0.908, 0.893}  //Ni+
};