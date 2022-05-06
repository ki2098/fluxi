#ifndef SOLVER_H
#define SOLVER_H

#include "setting.h"

void solver_sor(
    unsigned int  F[NNX][NNY][NNZ],
    double        P[NNX][NNY][NNZ],
    double       BP[NNX][NNY][NNZ][3],
    double       PS[NNX][NNY][NNZ],
    double        C[NNX][NNY][NNZ][6],
    double       &res
);

void solver_jacobi(
    unsigned int  F[NNX][NNY][NNZ],
    double        P[NNX][NNY][NNZ],
    double       PD[NNX][NNY][NNZ],
    double       BP[NNX][NNY][NNZ][3],
    double       PS[NNX][NNY][NNZ],
    double        C[NNX][NNY][NNZ][6],
    double       &res
);

#endif
