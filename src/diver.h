#ifndef DIVER_H
#define DIVER_H

#include "setting.h"

void diver(
    unsigned int   F[NNX][NNY][NNZ],
    double        UU[NNX][NNY][NNZ][3],
    double       DIV[NNX][NNY][NNZ],
    double         J[NNX][NNY][NNZ],
    double         &div
);

#endif
