#include "tp.h"

void _cell(
    unsigned int F[NNX][NNY][NNZ]
) {
//  set active cells
    for (int i = I0; i <= I1; i ++) {
        for (int j = J0; j <= J1; j ++) {
            for (int k = K0; k <= K1; k ++) {
                F[i][j][k] = f_set(F[i][j][k], ACTIVE, 1, MASK1);
            }
        }
    }

//  clear active cells in the square
    for (int i = I0sq; i <= I1sq; i ++) {
        for (int j = J0sq; j <= J1sq; j ++) {
            for (int k = K0sq; k <= K1sq; k ++) {
                F[i][j][k] = f_set(F[i][j][k], ACTIVE, 0, MASK1);
            }
        }
    }
}

void _face(
    unsigned int F[NNX][NNY][NNZ]
) {
    //  inflow and outflow boundary
    for (int j = J0; j <= J1; j ++) {
        for (int k = K0; k <= K1; k ++) {
            F[I0 - 1][j][k] = f_set(F[I0 - 1][j][k], F_E, INFLOW , MASK8);
            F[I1    ][j][k] = f_set(F[I1    ][j][k], F_E, OUTFLOW, MASK8);
        }
    }

//  slip boundary
    for (int i = I0; i <= I1; i ++) {
        for (int k = K0; k <= K1; k ++) {
            F[i][J0 - 1][k] = f_set(F[i][J0 - 1][k], F_N, SLIP, MASK8);
            F[i][J1    ][k] = f_set(F[i][J1    ][k], F_N, SLIP, MASK8);
        }
    }

//  square east and west boundary
    for (int j = J0sq; j <= J1sq; j ++) {
        for (int k = K0sq; k <= K1sq; k ++) {
            F[I0sq - 1][j][k] = f_set(F[I0sq - 1][j][k], F_E, WALL, MASK8);
            F[I1sq    ][j][k] = f_set(F[I1sq    ][j][k], F_E, WALL, MASK8);
        }
    }

//  square north and south boundary
    for (int i = I0sq; i <= I1sq; i ++) {
        for (int k = K0sq; k <= K1sq; k ++) {
            F[i][J0sq - 1][k] = f_set(F[i][J0sq - 1][k], F_N, WALL, MASK8);
            F[i][J1sq    ][k] = f_set(F[i][J1sq    ][k], F_N, WALL, MASK8);
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
                ac0 = f_see(F[i    ][j    ][k    ], ACTIVE, MASK1);
                ae1 = f_see(F[i + 1][j    ][k    ], ACTIVE, MASK1);
                an1 = f_see(F[i    ][j + 1][k    ], ACTIVE, MASK1);
                at1 = f_see(F[i    ][j    ][k + 1], ACTIVE, MASK1);
                if (!ac0) {
                    if (ae1) {
                        F[i][j][k] = f_set(F[i][j][k], M_E, 1, MASK1);
                    }
                    if (an1) {
                        F[i][j][k] = f_set(F[i][j][k], M_N, 1, MASK1);
                    }
                    if (at1) {
                        F[i][j][k] = f_set(F[i][j][k], M_T, 1, MASK1);
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
    double  G[NNX][NNY][NNZ][3],
    double  C[NNX][NNY][NNZ][6]
) {
    double dx, dy, dz;
    dx = (X1 - X0) / NX;
    dy = (Y1 - Y0) / NY;
    dz = (Z1 - Z0) / NZ;

    for (int i = 0; i < NNX; i ++) {
        for (int j = 0; j < NNY; j ++) {
            for (int k = 0; k < NNZ; k ++) {
                X[i][j][k][_X] = (i - I0) * dx + 0.5 * dx;
                X[i][j][k][_Y] = (j - J0) * dy + 0.5 * dy;
                X[i][j][k][_Z] = (k - K0) * dz + 0.5 * dz;
            }
        }
    }

    for (int i = I0 - 1; i <= I1 + 1; i ++) {
        for (int j = J0 - 1; j <= J1 + 1; j ++) {
            for (int k = K0 - 1; k <= K1 + 1; k ++) {
                double xc0, xe1, xw1;
                double yc0, yn1, ys1;
                double zc0, zt1, zb1;
                double x11, x22, x33;
                double xx1, xx2, xx3;
                double k11, k22, k33;
                double g11, g22, g33;
                double c01, c02, c03;
                double c07, c08, c09;
                double det;

                xc0 = X[i    ][j    ][k    ][_X];
                xe1 = X[i + 1][j    ][k    ][_X];
                xw1 = X[i - 1][j    ][k    ][_X];
                yc0 = X[i    ][j    ][k    ][_Y];
                yn1 = X[i    ][j + 1][k    ][_Y];
                ys1 = X[i    ][j - 1][k    ][_Y];
                zc0 = X[i    ][j    ][k    ][_Z];
                zt1 = X[i    ][j    ][k + 1][_Z];
                zb1 = X[i    ][j    ][k - 1][_Z];

                x11 = 0.5 * (xe1 - xw1);
                x22 = 0.5 * (yn1 - ys1);
                x33 = 0.5 * (zt1 - zb1);
                xx1 = xe1 - 2 * xc0 + xw1;
                xx2 = yn1 - 2 * yc0 + ys1;
                xx3 = zt1 - 2 * zc0 + zb1;
                k11 = 1 / x11;
                k22 = 1 / x22;
                k33 = 1 / x33;
                det = x11 * x22 * x33;
                c01 = 1 / (x11 * x11);
                c02 = 1 / (x22 * x22);
                c03 = 1 / (x33 * x33);
                c07 = xx1 / (x11 * x11 * x11);
                c08 = xx2 / (x22 * x22 * x22);
                c09 = xx3 / (x33 * x33 * x33);
                g11 = det * k11 * k11;
                g22 = det * k22 * k22;
                g33 = det * k33 * k33;

                KX[i][j][k][_X] = k11;
                KX[i][j][k][_Y] = k22;
                KX[i][j][k][_Z] = k33;
                C[ i][j][k][0 ] = c01;
                C[ i][j][k][1 ] = c02;
                C[ i][j][k][2 ] = c03;
                C[ i][j][k][3 ] = c07;
                C[ i][j][k][4 ] = c08;
                C[ i][j][k][5 ] = c09;
                G[ i][j][k][_X] = g11;
                G[ i][j][k][_Y] = g22;
                G[ i][j][k][_Z] = g33;
                J[ i][j][k]     = det;
            }
        }
    }
}