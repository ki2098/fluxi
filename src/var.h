#ifndef VAR_H
#define VAR_H

#include "setting.h"

/* flags */
unsigned int   F[NNX][NNY][NNZ]       = {0};

/* physical variables */
double         U[NNX][NNY][NNZ][3]    = {0};
double        UA[NNX][NNY][NNZ][3]    = {0};
double        UC[NNX][NNY][NNZ][3]    = {0};
double        UU[NNX][NNY][NNZ][3]    = {0};
double       UUA[NNX][NNY][NNZ][3]    = {0};
double        UT[NNX][NNY][NNZ][3]    = {0};
double         P[NNX][NNY][NNZ]       = {0};
double        PD[NNX][NNY][NNZ]       = {0};
double       DVR[NNX][NNY][NNZ]       = {0};
double       DVA[NNX][NNY][NNZ]       = {0};
double       SGS[NNX][NNY][NNZ]       = {0};
double        UR[NNX][NNY][NNZ][3]    = {0};
double        PR[NNX][NNY][NNZ]       = {0};

/* coordinate variables */
double         X[NNX][NNY][NNZ][3]    = {0};
double        KX[NNX][NNY][NNZ][3]    = {0};
double         J[NNX][NNY][NNZ]       = {0};
double         G[NNX][NNY][NNZ][3]    = {0};

double        HP[NNX][NNY][NNZ][3]    = {0};

/* boundary condition values */
double        BU[NNX][NNY][NNZ][3][3] = {0};
double        BP[NNX][NNY][NNZ][3]    = {0};

#endif
