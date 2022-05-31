#include "contra.h"
#include "bc.h"
#include "flag.h"

void contra(
    unsigned int  F[NNX][NNY][NNZ],
    double        U[NNX][NNY][NNZ][3],
    double       UC[NNX][NNY][NNZ][3],
    double       UU[NNX][NNY][NNZ][3],
    double       BU[NNX][NNY][NNZ][3][3],
    double        X[NNX][NNY][NNZ][3],
    double       KX[NNX][NNY][NNZ][3],
    double        J[NNX][NNY][NNZ]
) {
    #pragma acc kernels loop independent collapse(3) present(U, KX, UC, J)
    for (int i = I0 - 1; i <= I1 + 1; i ++) {
        for (int j = J0 - 1; j <= J1 + 1; j ++) {
            for (int k = K0 - 1; k <= K1 + 1; k ++) {
                double k1, k2, k3, de;
                double u, v, w;

                u   =  U[i][j][k][_U];
                v   =  U[i][j][k][_V];
                w   =  U[i][j][k][_W];
                k1  = KX[i][j][k][_X];
                k2  = KX[i][j][k][_Y];
                k3  = KX[i][j][k][_Z];
                de  =  J[i][j][k];

                UC[i][j][k][_U] = de * k1 * u;
                UC[i][j][k][_V] = de * k2 * v;
                UC[i][j][k][_W] = de * k3 * w;
            }
        }
    }

    #pragma acc kernels loop independent collapse(3) present(F, U, UC, UU, BU, X, KX, J) copyin(B)
    for (int i = I0 - 1; i <= I1; i ++) {
        for (int j = J0 - 1; j <= J1; j ++) {
            for (int k = K0 - 1; k <= K1; k ++) {
                unsigned int ff;
                unsigned int f3, m3, b3;
                double u0, u1, uf;
                double xk;
                double kx;
                double rf, di;
                ff = F[i][j][k];

                m3 = f_see(ff, _M_E, MASK1);
                f3 = f_see(ff, _F_E, MASK8);
                b3 = f_see(B[f3], _B_U, MASK2);
                if (f3) {
                    rf = U[i + m3][j][k][_U];
                    di = 0.5 - m3;
                    uf = bc_eva(b3, rf, di, BU[i][j][k][_E][_U]);
                    xk = X[i + 1][j][k][_X] -  X[i][j][k][_X];
                    kx = 1 / xk;
                    uf = 0.5 * (J[i][j][k] + J[i + 1][j][k]) * kx * uf;
                } else {
                    u0 = UC[i    ][j][k][_U];
                    u1 = UC[i + 1][j][k][_U];
                    uf = 0.5 * (u0 + u1); 
                }
                UU[i][j][k][_U] = uf;

                m3 = f_see(ff, _M_N, MASK1);
                f3 = f_see(ff, _F_N, MASK8);
                b3 = f_see(B[f3], _B_V, MASK2);
                if (f3) {
                    rf = U[i][j + m3][k][_V];
                    di = 0.5 - m3;
                    uf = bc_eva(b3, rf, di, BU[i][j][k][_N][_V]);
                    xk = X[i][j + 1][k][_Y] -  X[i][j][k][_Y];
                    kx = 1 / xk;
                    uf = 0.5 * (J[i][j][k] + J[i][j + 1][k]) * kx * uf;
                } else {
                    u0 = UC[i][j    ][k][_V];
                    u1 = UC[i][j + 1][k][_V];
                    uf = 0.5 * (u0 + u1); 
                }
                UU[i][j][k][_V] = uf;

                m3 = f_see(ff, _M_T, MASK1);
                f3 = f_see(ff, _F_T, MASK8);
                b3 = f_see(B[f3], _B_W, MASK2);
                if (f3) {
                    rf = U[i][j][k + m3][_W];
                    di = 0.5 - m3;
                    uf = bc_eva(b3, rf, di, BU[i][j][k][_T][_W]);
                    xk = X[i][j][k + 1][_Z] -  X[i][j][k][_Z];
                    kx = 1 / xk;
                    uf = 0.5 * (J[i][j][k] + J[i][j][k + 1]) * kx * uf;
                } else {
                    u0 = UC[i][j][k    ][_W];
                    u1 = UC[i][j][k + 1][_W];
                    uf = 0.5 * (u0 + u1);
                }
                UU[i][j][k][_W] = uf;
            }
        }
    }
}
