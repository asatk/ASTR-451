#define NSPECIES 10

enum species {HYDROGEN, CARBON, SODIUM, MAGNESIUM, SILICON, POTASSIUM, CALCIUM,
    CHROMIUM, IRON, NICKEL};

double G = 6.67e-11;    // gravitational constant, J m kg-2
double k = 1.38e-23;    // boltmann constant, J K-1
double mu = 1.66e-27;   // atomic mass unit, kg
double tsol = 5.777e3;  // effective temperature of the sun (Sol)
/**
 * ABUNDANCES (# atom/# hydrogen, in the Sun)
 * per Table 16.3 (Gray, 3ed)
 */
double A[] = {
    1.00e+0,    //H
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
double I[] = {
    13.598,     //H
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
double amu[] = {
    1.008,      //H
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
 * PARTITION FUNCTIONS for NEUTRAL SPECIES at theta=1.0
 * per Table D.2 (Gray, 3ed)
 */
double u0[] = {
    0.301,      //H
    0.964,      //C
    0.320,      //Na
    0.007,      //Mg
    0.976,      //Si
    0.351,      //K
    0.073,      //Ca
    1.022,      //Cr
    1.446,      //Fe
    1.490       //Ni
};

/**
 * PARTITION FUNCTIONS for SINGLY-IONIZIED SPECIES at theta=1.0
 * per Table D.2 (Gray, 3ed)
 */
double u1[] = {
    0.000,      //H+
    0.773,      //C+
    0.000,      //Na+
    0.301,      //Mg+
    0.755,      //Si+
    0.000,      //K+
    0.344,      //Ca+
    0.856,      //Cr+
    1.638,      //Fe+
    1.037       //Ni+
};