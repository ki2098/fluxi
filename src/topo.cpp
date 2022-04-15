#include "topo.h"
#include "topo2.h"
#include "boundary2.h"

void topo_init(
    int F[NNX][NNY][NNZ]
) {
//  set active cells
    for (int i = I0; i <= I1; i ++) {
        for (int j = J0; j <= J1; j ++) {
            for (int k = K0; k <= K1; k ++) {
                F[i][j][k] = set_bit(F[i][j][k], ACT, 1);
            }
        }
    }

//  clear active cells in the square
    for (int i = I0sq; i <= I1sq; i ++) {
        for (int j = J0sq; j <= J1sq; j ++) {
            for (int k = K0sq; k <= K1sq; k ++) {
                F[i][j][k] = set_bit(F[i][j][k], ACT, 0);
            }
        }
    }

//  inflow and outflow boundary
    for (int j = J0; j <= J1; j ++) {
        for (int k = K0; k <= K1; k ++) {
            F[I0 - 1][j][k] = set_face(F[I0 - 1][j][k], F_E, INFLOW);
            F[I1    ][j][k] = set_face(F[I1    ][j][k], F_E, OUTFLOW);
        }
    }

//  slip boundary
    for (int i = I0; i <= I1; i ++) {
        for (int k = K0; k <= K1; k ++) {
            F[i][J0 - 1][k] = set_face(F[i][J0 - 1][k], F_N, SLIP);
            F[i][J1    ][k] = set_face(F[i][J1    ][k], F_N, SLIP);
        }
    }

//  square east and west boundary
    for (int j = J0sq; j <= J1sq; j ++) {
        for (int k = K0sq; k <= K1sq; k ++) {
            F[I0sq - 1][j][k] = set_face(F[I0sq - 1][j][k], F_E, WALL);
            F[I1sq    ][j][k] = set_face(F[I1sq    ][j][k], F_E, WALL);
        }
    }

//  square north and south boundary
    for (int i = I0sq; i <= I1sq; i ++) {
        for (int k = K0sq; k <= K1sq; k ++) {
            F[i][J0sq - 1][k] = set_face(F[i][J0sq - 1][k], F_N, WALL);
            F[i][J1sq    ][k] = set_face(F[i][J1sq    ][k], F_N, WALL);
        }
    }
}