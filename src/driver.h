#ifndef DRIVER_H
#define DRIVER_H

#include "setting.h"

void driver_monitor(
    double UU[NNX][NNY][NNZ][3],
    double  X[NNX][NNY][NNZ][3],
    double KX[NNX][NNY][NNZ][3],
    double &avg
);

void driver_p_gradient(
    double  U[NNX][NNY][NNZ][3],
    double KX[NNX][NNY][NNZ][3],
    double  J[NNX][NNY][NNZ],
    double  ubar,
    double &driver_p
);

#endif
