#include <stdio.h>
#include <stdlib.h>
#include "tp.h"

void _cell(
    unsigned int F[NNX][NNY][NNZ]
) {
//  set active cells
    for (int i = I0; i <= I1; i ++) {
        for (int j = J0; j <= J1; j ++) {
            for (int k = K0; k <= K1; k ++) {
                F[i][j][k] = f_set(F[i][j][k], _ACTIVE, 1, MASK1);
            }
        }
    }

//  clear active cells in the square
    for (int i = I0sq; i <= I1sq; i ++) {
        for (int j = J0sq; j <= J1sq; j ++) {
            for (int k = K0sq; k <= K1sq; k ++) {
                F[i][j][k] = f_set(F[i][j][k], _ACTIVE, 0, MASK1);
            }
        }
    }
}

void _face(
    unsigned int F[NNX][NNY][NNZ]
) {
    // //  inflow and outflow boundary
    // for (int j = J0; j <= J1; j ++) {
    //     for (int k = K0; k <= K1; k ++) {
    //         F[I0 - 1][j][k] = f_set(F[I0 - 1][j][k], _F_E, INFLOW , MASK8);
    //         F[I1    ][j][k] = f_set(F[I1    ][j][k], _F_E, OUTFLOW, MASK8);
    //     }
    // }

    //  domain north and south boundary
    for (int i = I0; i <= I1; i ++) {
        for (int k = K0; k <= K1; k ++) {
            F[i][J0 - 1][k] = f_set(F[i][J0 - 1][k], _F_N, SLIP_Y, MASK8);
            F[i][J1    ][k] = f_set(F[i][J1    ][k], _F_N, SLIP_Y, MASK8);
        }
    }

    //  domain top and bottom boundary
    for (int i = I0; i <= I1; i ++) {
        for (int j = J0; j <= J1; j ++) {
            F[i][j][K0 - 1] = f_set(F[i][j][K0 - 1], _F_T, WALL, MASK8);
            F[i][j][K1    ] = f_set(F[i][j][K1    ], _F_T, WALL, MASK8);
        }
    }

    //  cube east and west boundary
    for (int j = J0sq; j <= J1sq; j ++) {
        for (int k = K0sq; k <= K1sq; k ++) {
            F[I0sq - 1][j][k] = f_set(F[I0sq - 1][j][k], _F_E, WALL, MASK8);
            F[I1sq    ][j][k] = f_set(F[I1sq    ][j][k], _F_E, WALL, MASK8);
        }
    }

    //  cube north and south boundary
    for (int i = I0sq; i <= I1sq; i ++) {
        for (int k = K0sq; k <= K1sq; k ++) {
            F[i][J0sq - 1][k] = f_set(F[i][J0sq - 1][k], _F_N, WALL, MASK8);
            F[i][J1sq    ][k] = f_set(F[i][J1sq    ][k], _F_N, WALL, MASK8);
        }
    }

    //  cube top boundary
    for (int i = I0sq; i <= I1sq; i ++) {
        for (int j = J0sq; j <= J1sq; j ++) {
            F[i][j][K1sq] = f_set(F[i][j][K1sq], _F_T, WALL, MASK8);
        }
    }
}

void _normal(
    unsigned int F[NNX][NNY][NNZ]
) {
    for (int i = I0 - 1; i <= I1; i ++) {
        for (int j = J0 - 1; j <= J1; j ++) {
            for (int k = K0 - 1; k <= K1; k ++) {
                int ac0, ae1, an1, at1;
                ac0 = f_see(F[i    ][j    ][k    ], _ACTIVE, MASK1);
                ae1 = f_see(F[i + 1][j    ][k    ], _ACTIVE, MASK1);
                an1 = f_see(F[i    ][j + 1][k    ], _ACTIVE, MASK1);
                at1 = f_see(F[i    ][j    ][k + 1], _ACTIVE, MASK1);
                if (!ac0) {
                    if (ae1 && f_see(F[i][j][k], _F_E, MASK8)) {
                        F[i][j][k] = f_set(F[i][j][k], _M_E, 1, MASK1);
                    }
                    if (an1 && f_see(F[i][j][k], _F_N, MASK8)) {
                        F[i][j][k] = f_set(F[i][j][k], _M_N, 1, MASK1);
                    }
                    if (at1 && f_see(F[i][j][k], _F_T, MASK8)) {
                        F[i][j][k] = f_set(F[i][j][k], _M_T, 1, MASK1);
                    }
                }
            }
        }
    }
}

void tp_f(
    unsigned int F[NNX][NNY][NNZ]
) {
    _cell(F);
    _face(F);
    _normal(F);
}

void tp_x(
    double  X[NNX][NNY][NNZ][3],
    double KX[NNX][NNY][NNZ][3],
    double  J[NNX][NNY][NNZ],
    double  G[NNX][NNY][NNZ][3]
) {
    FILE* f;
    char line[128];

    f = fopen("./mesh/x.cell", "r");
    for (int i = 0; i < NNX; i ++) {
        fgets(line, 128, f);
        double x = strtod(line, NULL);
        for (int j = 0; j < NNY; j ++) {
            for (int k = 0; k < NNZ; k ++) {
                X[i][j][k][_X] = x;
            }
        }
    }
    fclose(f);

    f = fopen("./mesh/y.cell", "r");
    for (int j = 0; j < NNY; j ++) {
        fgets(line, 128, f);
        double y = strtod(line, NULL);
        for (int i = 0; i < NNX; i ++) {
            for (int k = 0; k < NNZ; k ++) {
                X[i][j][k][_Y] = y;
            }
        }
    }
    fclose(f);

    f = fopen("./mesh/z.cell", "r");
    for (int k = 0; k < NNZ; k ++) {
        fgets(line, 128, f);
        double z = strtod(line, NULL);
        for (int i = 0; i < NNX; i ++) {
            for (int j = 0; j < NNY; j ++) {
                X[i][j][k][_Z] = z;
            }
        }
    }
    fclose(f);

    for (int i = I0 - 1; i <= I1 + 1; i ++) {
        for (int j = J0 - 1; j <= J1 + 1; j ++) {
            for (int k = K0 - 1; k <= K1 + 1; k ++) {
                double xe1, xw1;
                double yn1, ys1;
                double zt1, zb1;
                double x11, x22, x33;
                double k11, k22, k33;
                double g11, g22, g33;
                double det;

                xe1 = X[i + 1][j    ][k    ][_X];
                xw1 = X[i - 1][j    ][k    ][_X];
                yn1 = X[i    ][j + 1][k    ][_Y];
                ys1 = X[i    ][j - 1][k    ][_Y];
                zt1 = X[i    ][j    ][k + 1][_Z];
                zb1 = X[i    ][j    ][k - 1][_Z];

                x11 = 0.5 * (xe1 - xw1);
                x22 = 0.5 * (yn1 - ys1);
                x33 = 0.5 * (zt1 - zb1);
                k11 = 1 / x11;
                k22 = 1 / x22;
                k33 = 1 / x33;
                det = x11 * x22 * x33;
                g11 = det * k11 * k11;
                g22 = det * k22 * k22;
                g33 = det * k33 * k33;

                KX[i][j][k][_X] = k11;
                KX[i][j][k][_Y] = k22;
                KX[i][j][k][_Z] = k33;
                G[ i][j][k][_X] = g11;
                G[ i][j][k][_Y] = g22;
                G[ i][j][k][_Z] = g33;
                J[ i][j][k]     = det;
            }
        }
    }
}
