#include "ns.h"
#include "util.h"
#include <math.h>

double _ffvc_muscl(
    double uc0,
    double ue1,
    double ue2,
    double un1,
    double un2,
    double ut1,
    double ut2,
    double uw1,
    double uw2,
    double us1,
    double us2,
    double ub1,
    double ub2,
    double ufe,
    double vfn,
    double wft,
    double ufw,
    double vfs,
    double wfb,
    double det
) {
    double d1, d2, d3, d4;
    double s1, s2, s3, s4;
    double g1, g2, g3, g4, g5, g6;
    double u11, u10, u01, u00;
    double f1, f0;
    double adv;

    double kappa = 1.0 / 3.0;
    double beta  = (3.0 - kappa) / (1.0 - kappa);
    double m1    = 1.0 - kappa;
    double m2    = 1.0 + kappa;

    d4 = ue2 - ue1;
    d3 = ue1 - uc0;
    d2 = uc0 - uw1;
    d1 = uw1 - uw2;
    s4 = copysign(1.0, d4);
    s3 = copysign(1.0, d3);
    s2 = copysign(1.0, d2);
    s1 = copysign(1.0, d1);
    g6  = s4 * max(0.0, min(fabs(d4), s4 * beta * d3));
    g5  = s3 * max(0.0, min(fabs(d3), s3 * beta * d4));
    g4  = s3 * max(0.0, min(fabs(d3), s3 * beta * d2));
    g3  = s2 * max(0.0, min(fabs(d2), s2 * beta * d3));
    g2  = s2 * max(0.0, min(fabs(d2), s2 * beta * d1));
    g1  = s1 * max(0.0, min(fabs(d1), s1 * beta * d2));
    u11 = ue1 - 0.25 * (m1 * g6 + m2 * g5);
    u10 = uc0 + 0.25 * (m1 * g3 + m2 * g4);
    u01 = uc0 - 0.25 * (m1 * g4 + m2 * g3);
    u00 = uw1 + 0.25 * (m1 * g1 + m2 * g2);
    f1  = 0.5 * (ufe * (u11 + u10) - fabs(ufe) * (u11 - u10));
    f0  = 0.5 * (ufw * (u01 + u00) - fabs(ufw) * (u01 - u00));
    adv = (f1 - f0) / det;

    d4  = un2 - un1;
    d3  = un1 - uc0;
    d2  = uc0 - us1;
    d1  = us1 - us2;
    s4  = copysign(1.0, d4);
    s3  = copysign(1.0, d3);
    s2  = copysign(1.0, d2);
    s1  = copysign(1.0, d1);
    g6  = s4 * max(0.0, min(fabs(d4), s4 * beta * d3));
    g5  = s3 * max(0.0, min(fabs(d3), s3 * beta * d4));
    g4  = s3 * max(0.0, min(fabs(d3), s3 * beta * d2));
    g3  = s2 * max(0.0, min(fabs(d2), s2 * beta * d3));
    g2  = s2 * max(0.0, min(fabs(d2), s2 * beta * d1));
    g1  = s1 * max(0.0, min(fabs(d1), s1 * beta * d2));
    u11 = un1 - 0.25 * (m1 * g6 + m2 * g5);
    u10 = uc0 + 0.25 * (m1 * g3 + m2 * g4);
    u01 = uc0 - 0.25 * (m1 * g4 + m2 * g3);
    u00 = us1 + 0.25 * (m1 * g1 + m2 * g2);
    f1  = 0.5 * (vfn * (u11 + u10) - fabs(vfn) * (u11 - u10));
    f0  = 0.5 * (vfs * (u01 + u00) - fabs(vfs) * (u01 - u00));
    adv += (f1 - f0) / det;
}

void ns_correction_c(
    unsigned int   F[NNX][NNY][NNZ],
    double         U[NNX][NNY][NNZ][3],
    double        UD[NNX][NNY][NNZ][3],
    double        PP[NNX][NNY][NNZ],
    double        KX[NNX][NNY][NNZ][3],
    double         DT
) {
    for (int i = I0; i <= I1; i ++) {
        for (int j = J0; j <= J1; j ++) {
            for (int k = K0; k <= K1; k ++) {
                unsigned int fc0 = F[i][j][k];
                if (f_see(fc0, ACTIVE, MASK1)) {
                    unsigned int ffe, ffn, fft;
                    unsigned int ffw, ffs, ffb;
                    unsigned int mfe, mfn, mft;
                    unsigned int mfw, mfs, mfb;
                    unsigned int flw, fls, flb;

                    flw = F[i - 1][j    ][k    ];
                    fls = F[i    ][j - 1][k    ];
                    flb = F[i    ][j    ][k - 1];
                }
            }
        }
    }
}