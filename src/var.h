#ifndef VAR_H
#define VAR_H

#include "topo.h"

/* flags */
unsigned int   F[NNX][NNY][NNZ];

/* physical variables */
double         U[NNX][NNY][NNZ][3];
double        UA[NNX][NNY][NNZ][3];
double        UP[NNX][NNY][NNZ][3];
double        UD[NNX][NNY][NNZ][3];
double        UU[NNX][NNY][NNZ][3];
double       UUA[NNX][NNY][NNZ][3];
double       UUP[NNX][NNY][NNZ][3];
double       UUD[NNX][NNY][NNZ][3];
double         P[NNX][NNY][NNZ];
double        PD[NNX][NNY][NNZ];
double        PP[NNX][NNY][NNZ];
double       PPD[NNX][NNY][NNZ];
double       DVR[NNX][NNY][NNZ];
double       DVA[NNX][NNY][NNZ];
double       DVP[NNX][NNY][NNZ];
double       SGS[NNX][NNY][NNZ];

/* coordinate variables */
double         X[NNX][NNY][NNZ][3];
double        KX[NNX][NNY][NNZ][3];
double         J[NNX][NNY][NNZ];
double         G[NNX][NNY][NNZ][3];
double         C[NNX][NNY][NNZ][6];

/* boundary condition values */
double        BU[NNX][NNY][NNZ][3][3];
double        BP[NNX][NNY][NNZ][3];
double       BPP[NNX][NNY][NNZ][3];

/* runtime variables */
double         AD    = 0;
double         R     = 0;
int            ITER  = 0;
int            STEP  = 0;


#endif