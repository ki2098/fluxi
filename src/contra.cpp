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
    for (int i = I0 - 1; i <= I1 + 1; i ++) {
        for (int j = J0 - 1; j <= J1 + 1; j ++) {
            for (int k = K0 - 1; k <= K1 + 1; k ++) {
                double det, k1, k2, k3;
                double u, v, w;

                u   =  U[i][j][k][_U];
                v   =  U[i][j][k][_V];
                w   =  U[i][j][k][_W];
                k1  = KX[i][j][k][_X];
                k2  = KX[i][j][k][_Y];
                k3  = KX[i][j][k][_Z];
                det =  J[i][j][k];

                UC[i][j][k][_U] = det * k1 * u;
                UC[i][j][k][_V] = det * k2 * v;
                UC[i][j][k][_W] = det * k3 * w;
            }
        }
    }

    for (int i = I0 - 1; i <= I1; i ++) {
        for (int j = J0 - 1; j <= J1; j ++) {
            for (int k = K0 - 1; k <= K1; k ++) {
                unsigned int fla;
                unsigned int ffe, ffn, fft;
                unsigned int mfe, mfn, mft;
                double u0, u1, uf;
                double x1, x2, x3;
                double k1, k2, k3;
                double de;
                double rf, di, bc;
                unsigned int dr, nm;

                fla = F[i][j][k];
                ffe = f_see(fla, F_E, MASK8);
                ffn = f_see(fla, F_N, MASK8);
                fft = f_see(fla, F_T, MASK8);
                mfe = f_see(fla, M_E, MASK1);
                mfn = f_see(fla, M_N, MASK1);
                mft = f_see(fla, M_T, MASK1);

                if (ffe == NOTHING) {
                    u0 = UC[i    ][j][k][_U];
                    u1 = UC[i + 1][j][k][_U];
                    uf = 0.5 * (u0 + u1); 
                } else {
                    rf = U[i + mfe][j][k][_U];
                    di = 0.5 - mfe;
                    bc = BU[i][j][k][_X][_U];
                    dr = f_see(B[ffe], D_U, MASK1);
                    nm = f_see(B[ffe], N_U, MASK1);
                    uf = bc_evaluate(dr, nm, rf, di, bc);
                    x1 =         X[i + 1][j][k][_X] -  X[i][j][k][_X];
                    k2 = 0.5 * (KX[i + 1][j][k][_Y] + KX[i][j][k][_Y]);
                    k3 = 0.5 * (KX[i + 1][j][k][_Z] + KX[i][j][k][_Z]);
                    k1 = 1 / x1;
                    x2 = 1 / k2;
                    x3 = 1 / k3;
                    de = x1 * x2 * x3;
                    uf = de * k1 * uf;
                }
                UU[i][j][k][_U] = uf;

                if (ffn == NOTHING) {
                    u0 = UC[i][j    ][k][_V];
                    u1 = UC[i][j + 1][k][_V];
                    uf = 0.5 * (u0 + u1); 
                } else {
                    rf = U[i][j + mfn][k][_V];
                    di = 0.5 - mfn;
                    bc = BU[i][j][k][_Y][_V];
                    dr = f_see(B[ffn], D_V, MASK1);
                    nm = f_see(B[ffn], N_V, MASK1);
                    uf = bc_evaluate(dr, nm, rf, di, bc);
                    k1 = 0.5 * (KX[i][j + 1][k][_X] + KX[i][j][k][_X]);
                    x2 =         X[i][j + 1][k][_Y] -  X[i][j][k][_Y];
                    k3 = 0.5 * (KX[i][j + 1][k][_Z] + KX[i][j][k][_Z]);
                    x1 = 1 / k1;
                    k2 = 1 / x2;
                    x3 = 1 / k3;
                    de = x1 * x2 * x3;
                    uf = de * k2 * uf;
                }
                UU[i][j][k][_V] = uf;

                if (fft == NOTHING) {
                    u0 = UC[i][j][k    ][_W];
                    u1 = UC[i][j][k + 1][_W];
                    uf = 0.5 * (u0 + u1);
                } else {
                    rf = U[i][j][k + mft][_W];
                    di = 0.5 - mft;
                    bc = BU[i][j][k][_Z][_W];
                    dr = f_see(B[fft], D_W, MASK1);
                    nm = f_see(B[fft], N_W, MASK1);
                    uf = bc_evaluate(dr, nm, rf, di, bc);
                    k1 = 0.5 * (KX[i][j][k + 1][_X] + KX[i][j][k][_X]);
                    k2 = 0.5 * (KX[i][j][k + 1][_Y] + KX[i][j][k][_Y]);
                    x3 =         X[i][j][k + 1][_Z] -  X[i][j][k][_Z];
                    x1 = 1 / k1;
                    x2 = 1 / k2;
                    k3 = 1 / x3;
                    de = x1 * x2 * x3;
                    uf = de * k3 * uf;
                }
                UU[i][j][k][_W] = uf;
            }
        }
    }
}