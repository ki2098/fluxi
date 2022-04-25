#ifndef BC_H
#define BC_H

#include "setting.h"

#define B_PP 0

/* evaluate face value according to boundary condition */
static double bc_eva(unsigned int b, double reference, double distance, double value) {
    unsigned int drc = f_see(b, _BT_D, MASK1);
    unsigned int neu = f_see(b, _BT_N, MASK1);
    if (drc) {
        return value;
    }
    if (neu) {
        return reference + distance * value;
    }
    return 0;
}

/* pressure boundary value setting function */
void bc_u_init(
    double  U[NNX][NNY][NNZ][3],
    double BU[NNX][NNY][NNZ][3][3]
);

/* velocity boundary values setting function */
void bc_p_init(
    double  P[NNX][NNY][NNZ],
    double BP[NNX][NNY][NNZ][3]
);

/* slip velocity condition */
void bc_u_slip(
    double  U[NNX][NNY][NNZ][3],
    double BU[NNX][NNY][NNZ][3][3]
);

/* outflow velocity condition */
void bc_u_outflow(
    double  U[NNX][NNY][NNZ][3],
    double UU[NNX][NNY][NNZ][3],
    double BU[NNX][NNY][NNZ][3][3],
    double  X[NNX][NNY][NNZ][3],
    double KX[NNX][NNY][NNZ][3]
);

/* just copy it for 4 times */
void bc_u_periodic(
    double U[NNX][NNY][NNZ][3]
);

/* just copy it for 4 times */
void bc_p_periodic(
    double P[NNX][NNY][NNZ]
);

#endif
