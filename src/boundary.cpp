#include "boundary.h"
#include "boundary2.h"
#include "para.h"

/* slip velocity condition */
void _slip(
    double  U[NNX][NNY][NNZ][3],
    double BU[NNX][NNY][NNZ][3][3]
) {
    for (int i = I0; i <= I1; i ++) {
        for (int k = K0; k <= K1; k ++) {
            BU[i][J0 - 1][k][_N][_U] = U[i][J0][k][_U];
            BU[i][J0 - 1][k][_N][_V] = U[i][J0][k][_V];
            BU[i][J0 - 1][k][_N][_W] = U[i][J0][k][_W];
            BU[i][J1    ][k][_N][_U] = U[i][J1][k][_U];
            BU[i][J1    ][k][_N][_V] = U[i][J1][k][_V];
            BU[i][J1    ][k][_N][_W] = U[i][J1][k][_W];
        }
    }
}

/* outflow velocity condition */
void _outflow(
    double  U[NNX][NNY][NNZ][3],
    double UU[NNX][NNY][NNZ][3],
    double BU[NNX][NNY][NNZ][3][3],
    double  X[NNX][NNY][NNZ][3],
    double KX[NNX][NNY][NNZ][3],
    double  J[NNX][NNY][NNZ]
) {
    double UM = UINFLOW;
    for (int j = J0; j <= J1; j ++) {
        for (int k = K0; k <= K1; k ++) {
            double ufe, vfe, wfe;
            double ufw, vcc, wcc;
            double kfw, kfe, jfw;
            jfw = 0.5 * (J[I1    ][j][k]     + J[I1 - 1][j][k]);
            kfw =   1 / (X[I1    ][j][k][_X] - X[I1 - 1][j][k][_X]);
            kfe =   1 / (X[I1 + 1][j][k][_X] - X[I1    ][j][k][_X]);

            ufe = BU[I1    ][j][k][_X][_U];
            vfe = BU[I1    ][j][k][_X][_V];
            wfe = BU[I1    ][j][k][_X][_W];
            ufw = UU[I1 - 1][j][k][_U] / (jfw * kfw);
            vcc =  U[I1    ][j][k][_V];
            wcc =  U[I1    ][j][k][_W];

            ufe = ufe - DT * UM * kfe * (ufe - ufw);
            vfe = vfe - DT * UM * kfe * (vfe - vcc) * 2;
            wfe = wfe - DT * UM * kfe * (wfe - wcc) * 2;

            BU[I1    ][j][k][_X][_U] = ufe;
            BU[I1    ][j][k][_X][_V] = vfe;
            BU[I1    ][j][k][_X][_W] = wfe;
        }
    }
}

/* evaluate face value according to boundary condition */
double b_evaluate(
    unsigned int Drch, 
    unsigned int Neum, 
    double       reference, 
    double       distance, 
    double       value
) {
    if (Drch) {
        return value;
    }
    if (Neum) {
        return reference + distance * value;
    }
    return 0;
}

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
) {
    /* nothing to do */
}

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
) {
    if (timing == 1) {
        _outflow(U, UU, BU, X, KX, J);
    }
    else if (timing == 2) {
        _slip(U, BU);
    }
}

/* initialize boundary condition values*/
void b_init(
    double BU[NNX][NNY][NNZ][3][3],
    double BP[NNX][NNY][NNZ][3]
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