#include "bc.h"

/* what to be done at initialization */
void bc_u_init(
    double  U[NNX][NNY][NNZ][3],
    double BU[NNX][NNY][NNZ][3][3]
) {
    for (int j = J0; j <= J1; j ++) {
        for (int k = K0; k <= K1; k ++) {
            BU[I0 - 1][j][k][_E][_U] = UINFLOW;
            BU[I0 - 1][j][k][_E][_V] = VINFLOW;
            BU[I0 - 1][j][k][_E][_W] = WINFLOW;
            BU[I1    ][j][k][_E][_U] = UINFLOW;
            BU[I1    ][j][k][_E][_V] = VINFLOW;
            BU[I1    ][j][k][_E][_W] = WINFLOW;
        }
    }
}

/* what to be done at initialization */
void bc_p_init(
    double  P[NNX][NNY][NNZ],
    double BP[NNX][NNY][NNZ][3]
) {
    // nothing to do
}

/* outflow velocity condition */
void bc_u_outflow(
    double  U[NNX][NNY][NNZ][3],
    double UU[NNX][NNY][NNZ][3],
    double BU[NNX][NNY][NNZ][3][3],
    double  X[NNX][NNY][NNZ][3],
    double KX[NNX][NNY][NNZ][3]
) {
    double UM = UINFLOW;

    #pragma acc kernels loop independent collapse(2) present(U, UU, BU, X, KX) copyin(UM)
    for (int j = J0; j <= J1; j ++) {
        for (int k = K0; k <= K1; k ++) {
            double ufe, vfe, wfe;
            double ufw, vcc, wcc;
            double kfw, kfe, jfw;
            double kf1, kf2, kf3;
            double xf1, xf2, xf3;

            xf1 =         X[I1][j][k][_X] -  X[I1 - 1][j][k][_X];
            kf2 = 0.5 * (KX[I1][j][k][_Y] + KX[I1 - 1][j][k][_Y]);
            kf3 = 0.5 * (KX[I1][j][k][_Z] + KX[I1 - 1][j][k][_Z]);
            kf1 = 1 / xf1;
            xf2 = 1 / kf2;
            xf3 = 1 / kf3;
            jfw = xf1 * xf2 * xf3;

            kfw = kf1;
            kfe = 1 / (X[I1 + 1][j][k][_X] - X[I1][j][k][_X]);

            ufe = BU[I1    ][j][k][_E][_U];
            vfe = BU[I1    ][j][k][_E][_V];
            wfe = BU[I1    ][j][k][_E][_W];
            ufw = UU[I1 - 1][j][k][_U] / (jfw * kfw);
            vcc =  U[I1    ][j][k][_V];
            wcc =  U[I1    ][j][k][_W];

            ufe = ufe - DT * UM * kfe * (ufe - ufw);
            vfe = vfe - DT * UM * kfe * (vfe - vcc) * 2;
            wfe = wfe - DT * UM * kfe * (wfe - wcc) * 2;

            BU[I1][j][k][_E][_U] = ufe;
            BU[I1][j][k][_E][_V] = vfe;
            BU[I1][j][k][_E][_W] = wfe;
        }
    }
}
