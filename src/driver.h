#ifndef DRIVER_H
#define DRIVER_H

#include "setting.h"

void driver_monitor(
    double  U[NNX][NNY][NNZ][3],
    double  J[NNX][NNY][NNZ],
    double &avg
);

void driver_p_gradient(
    double  U[NNX][NNY][NNZ][3],
    double  J[NNX][NNY][NNZ],
    double  ubar,
    double &driver_p
);

#endif
