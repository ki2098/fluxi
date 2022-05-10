#include "bc.h"

/* what to be done at initialization */
void bc_u_init(
    double  U[NNX][NNY][NNZ][3],
    double BU[NNX][NNY][NNZ][3][3]
) {
    // for (int j = J0; j <= J1; j ++) {
    //     for (int k = K0; k <= K1; k ++) {
    //         BU[I0 - 1][j][k][_E][_U] = UINFLOW;
    //         BU[I0 - 1][j][k][_E][_V] = VINFLOW;
    //         BU[I0 - 1][j][k][_E][_W] = WINFLOW;
    //         BU[I1    ][j][k][_E][_U] = UINFLOW;
    //         BU[I1    ][j][k][_E][_V] = VINFLOW;
    //         BU[I1    ][j][k][_E][_W] = WINFLOW;
    //     }
    // }

    // nothing to do
}

/* what to be done at initialization */
void bc_p_init(
    double  P[NNX][NNY][NNZ],
    double BP[NNX][NNY][NNZ][3]
) {
    // nothing to do
}

// /* outflow velocity condition */
// void bc_u_outflow(
//     double  U[NNX][NNY][NNZ][3],
//     double UU[NNX][NNY][NNZ][3],
//     double BU[NNX][NNY][NNZ][3][3],
//     double  X[NNX][NNY][NNZ][3],
//     double KX[NNX][NNY][NNZ][3]
// ) {
//     double UM = UINFLOW;

//     #pragma acc kernels loop independent collapse(2) present(U, UU, BU, X, KX) copyin(UM)
//     for (int j = J0; j <= J1; j ++) {
//         for (int k = K0; k <= K1; k ++) {
//             double ufe, vfe, wfe;
//             double ufw, vcc, wcc;
//             double kfw, kfe, jfw;
//             double kf1, kf2, kf3;
//             double xf1, xf2, xf3;

//             xf1 =         X[I1][j][k][_X] -  X[I1 - 1][j][k][_X];
//             kf2 = 0.5 * (KX[I1][j][k][_Y] + KX[I1 - 1][j][k][_Y]);
//             kf3 = 0.5 * (KX[I1][j][k][_Z] + KX[I1 - 1][j][k][_Z]);
//             kf1 = 1 / xf1;
//             xf2 = 1 / kf2;
//             xf3 = 1 / kf3;
//             jfw = xf1 * xf2 * xf3;

//             kfw = kf1;
//             kfe = 1 / (X[I1 + 1][j][k][_X] - X[I1][j][k][_X]);

//             ufe = BU[I1    ][j][k][_E][_U];
//             vfe = BU[I1    ][j][k][_E][_V];
//             wfe = BU[I1    ][j][k][_E][_W];
//             ufw = UU[I1 - 1][j][k][_U] / (jfw * kfw);
//             vcc =  U[I1    ][j][k][_V];
//             wcc =  U[I1    ][j][k][_W];

//             ufe = ufe - DT * UM * kfe * (ufe - ufw);
//             vfe = vfe - DT * UM * kfe * (vfe - vcc) * 2;
//             wfe = wfe - DT * UM * kfe * (wfe - wcc) * 2;

//             BU[I1][j][k][_E][_U] = ufe;
//             BU[I1][j][k][_E][_V] = vfe;
//             BU[I1][j][k][_E][_W] = wfe;
//         }
//     }
// }

void bc_u_periodic(
    double  U[NNX][NNY][NNZ][3]
) {
    #pragma acc kernels loop independent collapse(2) present(U)
    for (int j = J0; j <= J1; j ++) {
        for (int k = K0; k <= K1; k ++) {
            U[I0 - 1][j][k][_U] = U[I1    ][j][k][_U];
            U[I0 - 1][j][k][_V] = U[I1    ][j][k][_V];
            U[I0 - 1][j][k][_W] = U[I1    ][j][k][_W];
            U[I0 - 2][j][k][_U] = U[I1 - 1][j][k][_U];
            U[I0 - 2][j][k][_V] = U[I1 - 1][j][k][_V];
            U[I0 - 2][j][k][_W] = U[I1 - 1][j][k][_W];
            U[I1 + 1][j][k][_U] = U[I0    ][j][k][_U];
            U[I1 + 1][j][k][_V] = U[I0    ][j][k][_V];
            U[I1 + 1][j][k][_W] = U[I0    ][j][k][_W];
            U[I1 + 2][j][k][_U] = U[I0 + 1][j][k][_U];
            U[I1 + 2][j][k][_V] = U[I0 + 1][j][k][_V];
            U[I1 + 2][j][k][_W] = U[I0 + 1][j][k][_W];
        }
    }

    #pragma acc kernels loop independent collapse(2) present(U)
    for (int i = I0; i <= I1; i ++) {
        for (int k = K0; k <= K1; k ++) {
            U[i][J0 - 1][k][_U] = U[i][J1    ][k][_U];
            U[i][J0 - 1][k][_V] = U[i][J1    ][k][_V];
            U[i][J0 - 1][k][_W] = U[i][J1    ][k][_W];
            U[i][J0 - 2][k][_U] = U[i][J1 - 1][k][_U];
            U[i][J0 - 2][k][_V] = U[i][J1 - 1][k][_V];
            U[i][J0 - 2][k][_W] = U[i][J1 - 1][k][_W];
            U[i][J1 + 1][k][_U] = U[i][J0    ][k][_U];
            U[i][J1 + 1][k][_V] = U[i][J0    ][k][_V];
            U[i][J1 + 1][k][_W] = U[i][J0    ][k][_W];
            U[i][J1 + 2][k][_U] = U[i][J0 + 1][k][_U];
            U[i][J1 + 2][k][_V] = U[i][J0 + 1][k][_V];
            U[i][J1 + 2][k][_W] = U[i][J0 + 1][k][_W];
        }
    }
}

void bc_p_driver(
    double P[NNX][NNY][NNZ],
    double driver_p
) {
    #pragma acc kernels loop independent collapse(2) present(P) copyin(driver_p)
    for (int j = J0; j <= J1; j ++) {
        for (int k = K0; k <= K1; k ++) {
            P[I0 - 1][j][k] = P[I1][j][k] + driver_p;
            P[I1 + 1][j][k] = P[I0][j][k] - driver_p;
        }
    }
}

void bc_p_periodic(
    double P[NNX][NNY][NNZ]
) {
    #pragma acc kernels loop independent collapse(2) present(P)
    for (int i = I0; i <= I1; i ++) {
        for (int k = K0; k <= K1; k ++) {
            P[i][J0 - 1][k] = P[i][J1][k];
            P[i][J1 + 1][k] = P[i][J0][k];
        }
    }
}
