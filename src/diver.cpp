#include <math.h>
#include "diver.h"

void diver(
    unsigned int   F[NNX][NNY][NNZ],
    double        UU[NNX][NNY][NNZ][3],
    double       DIV[NNX][NNY][NNZ],
    double         J[NNX][NNY][NNZ],
    double         &div
) {
    double sum = 0;
    int    cnt = 0;

    #pragma acc kernels loop independent collapse(3) reduction(+:sum, cnt) present(F, UU, DIV, J) copy(sum, cnt)
    for (int i = I0; i <= I1; i ++) {
        for (int j = J0; j <= J1; j ++) {
            for (int k = K0; k <= K1; k ++) {
                if (f_see(F[i][j][k], _ACTIVE, MASK1)) {
                    double ufe, vfn, wft;
                    double ufw, vfs, wfb;
                    double det;
                    double dc0;

                    ufe = UU[i    ][j    ][k    ][_U];
                    vfn = UU[i    ][j    ][k    ][_V];
                    wft = UU[i    ][j    ][k    ][_W];
                    ufw = UU[i - 1][j    ][k    ][_U];
                    vfs = UU[i    ][j - 1][k    ][_V];
                    wfb = UU[i    ][j    ][k - 1][_W];
                    det =  J[i    ][j    ][k    ];

                    dc0 = (ufe + vfn + wft - ufw - vfs - wfb) / det;
                    sum = sum + dc0 * dc0;
                
                    DIV[i][j][k] = dc0;
                    cnt += 1;
                }
                
            }
        }
    }

    div = sqrt(sum / cnt);
}
