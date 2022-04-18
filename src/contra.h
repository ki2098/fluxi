#ifndef CONTRA_H
#define CONTRA_H

#include "setting.h"

void contra(
    unsigned int  F[NNX][NNY][NNZ],
    double        U[NNX][NNY][NNZ][3],
    double       UC[NNX][NNY][NNZ][3],
    double       UU[NNX][NNY][NNZ][3],
    double       BU[NNX][NNY][NNZ][3][3],
    double        X[NNX][NNY][NNZ][3],
    double       KX[NNX][NNY][NNZ][3],
    double        J[NNX][NNY][NNZ]
)

#endif