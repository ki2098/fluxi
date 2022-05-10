#include "driver.h"

void driver_monitor(
    double UU[NNX][NNY][NNZ][3],
    double  X[NNX][NNY][NNZ][3],
    double KX[NNX][NNY][NNZ][3],
    double &avg
) {
    double sum = 0;
    int    cnt = 0;
    #pragma acc kernels loop independent collapse(2) reduction(+:sum, cnt) present(UU, X, KX) copy(sum, cnt)
    for (int j = J0; j <= J1; j ++) {
        for (int k = K0 ; k <= K1; k ++) {
            double uu, x1, x2, x3, k1, k2, k3, de;

            uu = UU[I0][j][k][_U];
            x1 =         X[I0][j][k][_X] -  X[I0 - 1][j][k][_X];
            k2 = 0.5 * (KX[I0][j][k][_Y] + KX[I0 - 1][j][k][_Y]);
            k3 = 0.5 * (KX[I0][j][k][_Z] + KX[I0 - 1][j][k][_Z]);
            k1 = 1 / x1;
            x2 = 1 / k2;
            x3 = 1 / k3;
            de = x1 * x2 * x3;

            uu = uu / (de * k1);
            
            cnt += 1;
            sum += uu;
        }
    }
    avg = sum / cnt;
}

void driver_p_gradient(
    double  U[NNX][NNY][NNZ][3],
    double KX[NNX][NNY][NNZ][3],
    double  J[NNX][NNY][NNZ],
    double  ubar,
    double &driver_p
) {
    double m0 = 0;
    double m1 = 0;
    double a  = 0;
    #pragma acc kernels loop independent collapse(2) reduction(+:m0, m1, a) present(U, J, KX) copyin(ubar) copy(m0, m1, a)
    for (int i = I0; i <= I1; i ++) {
        for (int j = J0; j <= J1; j ++) {
            for (int k = K0; k <= K1; k ++) {
                m0 += U[i][j][k][_U] * J[i][j][k];
                m1 += ubar           * J[i][j][k];
                if (i == I0) {
                    a += 1.0 / (KX[i][j][k][_Y] * KX[i][j][k][_Z]);
                }
            }
        }
    }
    driver_p = (m1 - m0) / (DT * a);
}
