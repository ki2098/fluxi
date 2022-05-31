#include "driver.h"

void driver_monitor(
    double  U[NNX][NNY][NNZ][3],
    double  J[NNX][NNY][NNZ],
    double &avg
) {
    double m0 = 0;
    double m1 = 0;
    #pragma acc kernels loop independent collapse(2) reduction(+:m0, m1) present(U, J) copy(m0, m1)
    for (int j = J0; j <= J1; j ++) {
        for (int k = K0; k <= K1; k ++) {
            m0 += U[I0][j][k][_U] * J[I0][j][k];
            m1 += J[I0][j][k];
        }
    }
    avg = m0 / m1;
}

void driver_p_gradient(
    double  U[NNX][NNY][NNZ][3],
    double  J[NNX][NNY][NNZ],
    double  ubar,
    double &driver_p
) {
    double m0 = 0;
    double m1 = 0;
    #pragma acc kernels loop independent collapse(2) reduction(+:m0, m1) present(U, J) copyin(ubar) copy(m0, m1)
    for (int i = I0; i <= I1; i ++) {
        for (int j = J0; j <= J1; j ++) {
            for (int k = K0; k <= K1; k ++) {
                m0 += U[i][j][k][_U] * J[i][j][k];
                m1 += ubar           * J[i][j][k];
            }
        }
    }
    double a = (Z1 - Z0) * (Y1 - Y0);
    driver_p = (m1 - m0) / (DT * a);
}
