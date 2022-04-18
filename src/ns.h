#ifndef NS_H
#define NS_H

#include "setting.h"

void ns_pseudo_c(
    unsigned int   F[NNX][NNY][NNZ],
    double         U[NNX][NNY][NNZ][3],
    double        UA[NNX][NNY][NNZ][3],
    double        UU[NNX][NNY][NNZ][3],
    double        BU[NNX][NNY][NNZ][3][3],
    double       SGS[NNX][NNY][NNZ],
    double        KX[NNX][NNY][NNZ][3],
    double         J[NNX][NNY][NNZ],
    double         C[NNX][NNY][NNZ][6],
    double         ALPHA,
    double         DT,
    double         RI
);

void ns_correction_c(
    unsigned int   F[NNX][NNY][NNZ],
    double         U[NNX][NNY][NNZ][3],
    double        UD[NNX][NNY][NNZ][3],
    double        BU[NNX][NNY][NNZ][3][3],
    double        PP[NNX][NNY][NNZ],
    double        KX[NNX][NNY][NNZ][3],
    double         DT
);

void ns_correction_f(
    unsigned int   F[NNX][NNY][NNZ],
    double         U[NNX][NNY][NNZ][3],
    double        UU[NNX][NNY][NNZ][3],
    double       UUD[NNX][NNY][NNZ][3],
    double        BU[NNX][NNY][NNZ][3][3],
    double        PP[NNX][NNY][NNZ],
    double        KX[NNX][NNY][NNZ][3],
    double         G[NNX][NNY][NNZ],
    double         DT
);

#endif