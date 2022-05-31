#ifndef SMAC_H
#define SMAC_H

#include "setting.h"

void ns_pseudo_c(
    unsigned int   F[NNX][NNY][NNZ],
    double         U[NNX][NNY][NNZ][3],
    double        UA[NNX][NNY][NNZ][3],
    double        UU[NNX][NNY][NNZ][3],
    double        UT[NNX][NNY][NNZ][3],
    double        BU[NNX][NNY][NNZ][3][3],
    double       SGS[NNX][NNY][NNZ],
    double         X[NNX][NNY][NNZ][3],
    double         J[NNX][NNY][NNZ],
    double         G[NNX][NNY][NNZ][3]
);

void ns_correction_c(
    unsigned int   F[NNX][NNY][NNZ],
    double         U[NNX][NNY][NNZ][3],
    double        UD[NNX][NNY][NNZ][3],
    double         P[NNX][NNY][NNZ],
    double        BP[NNX][NNY][NNZ][3],
    double        KX[NNX][NNY][NNZ][3]
);

void ns_correction_f(
    unsigned int   F[NNX][NNY][NNZ],
    double         U[NNX][NNY][NNZ][3],
    double        UU[NNX][NNY][NNZ][3],
    double       UUD[NNX][NNY][NNZ][3],
    double        BU[NNX][NNY][NNZ][3][3],
    double         P[NNX][NNY][NNZ],
    double        BP[NNX][NNY][NNZ][3],
    double         X[NNX][NNY][NNZ][3],
    double        KX[NNX][NNY][NNZ][3],
    double         G[NNX][NNY][NNZ][3]
);

void ns_correction_p(
    double  P[NNX][NNY][NNZ],
    double PP[NNX][NNY][NNZ]
);

#endif
