#ifndef TURB_H
#define TURB_H

#include "setting.h"

void turb_smagorinsky(
    unsigned int   F[NNX][NNY][NNZ],
    double         U[NNX][NNY][NNZ][3],
    double        BU[NNX][NNY][NNZ][3][3],
    double         X[NNX][NNY][NNZ][3],
    double        KX[NNX][NNY][NNZ][3],
    double         J[NNX][NNY][NNZ],
    double       SGS[NNX][NNY][NNZ]
);

#endif
