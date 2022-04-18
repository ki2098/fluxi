#include "bc.h"

/* what to be done at initialization */
void _bc_u_init(
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
    for (int i = I0; i <= I1; i ++) {
        for (int k = K0; k <= K1; k ++) {
            BU[i][J0 - 1][k][_N][_U] = UINIT;
            BU[i][J0 - 1][k][_N][_V] = 0;
            BU[i][J0 - 1][k][_N][_W] = 0;
            BU[i][J1    ][k][_N][_U] = UINIT;
            BU[i][J1    ][k][_N][_V] = 0;
            BU[i][J1    ][k][_N][_W] = 0;
        }
    }
    for (int i = 0; i < NNX; i ++) {
        for (int j = 0; j < NNY; j ++) {
            U[i][j][K0 - 1][_U] = U[i][j][K0][_U];
            U[i][j][K0 - 1][_V] = U[i][j][K0][_V];
            U[i][j][K0 - 1][_W] = U[i][j][K0][_W];
            U[i][j][K0 - 2][_U] = U[i][j][K0][_U];
            U[i][j][K0 - 2][_V] = U[i][j][K0][_V];
            U[i][j][K0 - 2][_W] = U[i][j][K0][_W];
            U[i][j][K1 + 1][_U] = U[i][j][K1][_U];
            U[i][j][K1 + 1][_V] = U[i][j][K1][_V];
            U[i][j][K1 + 1][_W] = U[i][j][K1][_W];
            U[i][j][K1 + 2][_U] = U[i][j][K1][_U];
            U[i][j][K1 + 2][_V] = U[i][j][K1][_V];
            U[i][j][K1 + 2][_W] = U[i][j][K1][_W];
        }
    }
}

/* what to be done at initialization */
void _bc_p_init(
    double  P[NNX][NNY][NNZ],
    double BP[NNX][NNY][NNZ][3]
) {
    for (int i = 0; i < NNX; i ++) {
        for (int j = 0; j < NNY; j ++) {
            P[i][j][K0 - 1] = P[i][j][K0];
            P[i][j][K0 - 2] = P[i][j][K0];
            P[i][j][K1 + 1] = P[i][j][K1];
            P[i][j][K1 + 2] = P[i][j][K1];
        }
    }
}

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

/* just copy it for 4 times */
void _periodic_u(
    double U[NNX][NNY][NNZ][3]
) {
    for (int i = 0; i < NNX; i ++) {
        for (int j = 0; j < NNY; j ++) {
            U[i][j][K0 - 1][_U] = U[i][j][K0][_U];
            U[i][j][K0 - 1][_V] = U[i][j][K0][_V];
            U[i][j][K0 - 1][_W] = U[i][j][K0][_W];
            U[i][j][K0 - 2][_U] = U[i][j][K0][_U];
            U[i][j][K0 - 2][_V] = U[i][j][K0][_V];
            U[i][j][K0 - 2][_W] = U[i][j][K0][_W];
            U[i][j][K1 + 1][_U] = U[i][j][K1][_U];
            U[i][j][K1 + 1][_V] = U[i][j][K1][_V];
            U[i][j][K1 + 1][_W] = U[i][j][K1][_W];
            U[i][j][K1 + 2][_U] = U[i][j][K1][_U];
            U[i][j][K1 + 2][_V] = U[i][j][K1][_V];
            U[i][j][K1 + 2][_W] = U[i][j][K1][_W];
        }
    }
}

/* just copy it for 4 times */
void _periodic_p(
    double P[NNX][NNY][NNZ]
) {
    for (int i = 0; i < NNX; i ++) {
        for (int j = 0; j < NNY; j ++) {
            P[i][j][K0 - 1] = P[i][j][K0];
            P[i][j][K0 - 2] = P[i][j][K0];
            P[i][j][K1 + 1] = P[i][j][K1];
            P[i][j][K1 + 2] = P[i][j][K1];
        }
    }
}

/* pressure boundary value setting function */
void bc_p(
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
/*  
    timing: 
    0 if for initialization
    2 is for after each iteration in the linear solver
    3 is reserved for periodic conditions
*/
/*
    π is only subject to peridic conditions
    since all conditions in p lead to 0 value or 0 gradient of π
    it doesn't need to be initialized either
    because at the beginning, it is surely all 0
*/
    if (timing == 0) {
        _bc_p_init(P, BP);
    } else if (timing == 2) {
        /* nothing to do */
    } else if (timing == 3) {
        _periodic_p(P);
    }
}

/* velocity boundary values setting function */
void bc_u(
    double  U[NNX][NNY][NNZ][3],
    double UU[NNX][NNY][NNZ][3],
    double BU[NNX][NNY][NNZ][3][3],
    double  X[NNX][NNY][NNZ][3],
    double KX[NNX][NNY][NNZ][3],
    double  J[NNX][NNY][NNZ],
    double  C[NNX][NNY][NNZ][6],
    int     timing
) {
/*  
    timing: 
    0 is for initialization
    1 is for just after u*@center is calculated
    2 is for after projection from uP@center to u@center
    3 is reserved for periodic conditions
*/
    if (timing == 0) {
        _bc_u_init(U, BU);
    } else if (timing == 1) {
        _outflow(U, UU, BU, X, KX, J);
    } else if (timing == 2) {
        /* nothing to do */
    } else if (timing == 3) {
        _periodic_u(U);
    }
}
