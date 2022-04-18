#ifndef TP_H
#define TP_H

#include "setting.h"

void tp_f(
    unsigned int F[NNX][NNY][NNZ]
);

void tp_x(
    double  X[NNX][NNY][NNZ][3],
    double KX[NNX][NNY][NNZ][3],
    double  J[NNX][NNY][NNZ],
    double  G[NNX][NNY][NNZ][3],
    double  C[NNX][NNY][NNZ][6]
);

#endif