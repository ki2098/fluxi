#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "topo.h"

/* evaluate face value according to boundary condition */
double b_evaluate(
    unsigned int Drch, 
    unsigned int Neum, 
    double       reference, 
    double       distance, 
    double       value
);

/* pressure boundary value setting function */
void b_p(
    double  P[NNX][NNY][NNZ],
    double BP[NNX][NNY][NNZ][3],
    double  U[NNX][NNY][NNZ][3],
    double UU[NNX][NNY][NNZ][3],
    double  X[NNX][NNY][NNZ][3],
    double KX[NNX][NNY][NNZ][3],
    double  J[NNX][NNY][NNZ],
    double  C[NNX][NNY][NNZ][6],
    int     timing
);

/* velocity boundary values setting function */
void b_u(
    double  U[NNX][NNY][NNZ][3],
    double UU[NNX][NNY][NNZ][3],
    double BU[NNX][NNY][NNZ][3][3],
    double  X[NNX][NNY][NNZ][3],
    double KX[NNX][NNY][NNZ][3],
    double  J[NNX][NNY][NNZ],
    double  C[NNX][NNY][NNZ][6],
    int     timing
);

/* initialize boundary condition values*/
void b_init(
    double BU[NNX][NNY][NNZ][3][3],
    double BP[NNX][NNY][NNZ][3]
)

#endif