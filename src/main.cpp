#include <stdio.h>
#include "tp.h"
#include "bc.h"
#include "var.h"

/* runtime variables */
double         AD    = 0;
double         R     = 0;
int            ITER  = 0;
int            STEP  = 0;

void _param_out(void) {
    FILE *fo;
    fo = fopen("para.csv", "w+t");
    if (fo == NULL) {
        printf("\nERROR when opening file\n");
        fflush(stdout);
    }
    else {
        fprintf(fo, "x,y,z,active,fe,fn,ft,me,mn,mt,k1x1,k2x2,k3x3,j,c1,c2,c3,c7,c8,c9,g11,g22,g33\n");
        double x, y, z;
        int    flag, active;
        int    fe, fn, ft;
        int    me, mn, mt;
        double k1x1, k2x2, k3x3;
        double det;
        double c1, c2, c3, c7, c8, c9;
        double g11, g22, g33;
        for (int k = 0; k < NNZ; k ++) {
            for (int j = 0; j < NNY; j ++) {
                for (int i = 0; i < NNX; i ++) {
                    x      =  X[i][j][k][_X];
                    y      =  X[i][j][k][_Y];
                    z      =  X[i][j][k][_Z];
                    flag   =  F[i][j][k];
                    k1x1   = KX[i][j][k][_X];
                    k2x2   = KX[i][j][k][_Y];
                    k3x3   = KX[i][j][k][_Z];
                    det    =  J[i][j][k];
                    c1     =  C[i][j][k][0 ];
                    c2     =  C[i][j][k][1 ];
                    c3     =  C[i][j][k][2 ];
                    c7     =  C[i][j][k][3 ];
                    c8     =  C[i][j][k][4 ];
                    c9     =  C[i][j][k][5 ];
                    g11    =  G[i][j][k][0 ];
                    g22    =  G[i][j][k][1 ];
                    g33    =  G[i][j][k][2 ];
                    active = f_see(flag, ACTIVE, MASK1);
                    fe     = f_see(flag, F_E   , MASK8);
                    fn     = f_see(flag, F_N   , MASK8);
                    ft     = f_see(flag, F_T   , MASK8);
                    me     = f_see(flag, M_E   , MASK1);
                    mn     = f_see(flag, M_N   , MASK1);
                    mt     = f_see(flag, M_T   , MASK1);
                    fprintf(fo, "%lf,%lf,%lf,%d,%d,%d,%d,%d,%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", x, y, z, active, fe, fn, ft, me, mn, mt, k1x1, k2x2, k3x3, det, c1, c2, c3, c7, c8, c9, g11, g22, g33);
                }
            }
        }
    }
}

void _var_out(void) {
    FILE *fo;
    fo = fopen("var.csv", "w+t");
    if (fo == NULL) {
        printf("\nERROR when opening file\n");
        fflush(stdout);
    }
    else {
        fprintf(fo, "x,y,z,u,v,w,p,dvr\n");
        double x, y, z, u, v, w, p, dvr;
        for (int k = K0; k <= K1; k ++) {
            for (int j = J0; j <= J1; j ++) {
                for (int i = I0; i <= I1; i ++) {
                    x   =   X[i][j][k][_X];
                    y   =   X[i][j][k][_Y];
                    z   =   X[i][j][k][_Z];
                    u   =   U[i][j][k][_U];
                    v   =   U[i][j][k][_V];
                    w   =   U[i][j][k][_W];
                    p   =   P[i][j][k];
                    dvr = DVR[i][j][k];
                    fprintf(fo, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", x, y, z, u, v, w, p, dvr);
                }
            }
        }
    }
}

void _boundary_out(void) {
    printf("nu du nv dv nw dw np dp\n");
    for (int i = 0; i <= NB; i ++) {
        unsigned int b  = B[i];
        unsigned int nu = f_see(b, N_U, MASK1);
        unsigned int du = f_see(b, D_U, MASK1);
        unsigned int nv = f_see(b, N_V, MASK1);
        unsigned int dv = f_see(b, D_V, MASK1);
        unsigned int nw = f_see(b, N_W, MASK1);
        unsigned int dw = f_see(b, D_W, MASK1);
        unsigned int np = f_see(b, N_P, MASK1);
        unsigned int dp = f_see(b, D_P, MASK1);
        printf("%u  %u  %u  %u  %u  %u  %u  %u\n", nu, du, nv, dv, nw, dw, np, dp);
    }
}

void _init(void) {
    tp_x(X, KX, J, G, C);
    tp_f(F);
    
    for (int i = I0; i <= I1; i ++) {
        for (int j = J0; j <= J1; j ++) {
            for (int k = K0; k <= K1; k ++) {
                if (f_see(F[i][j][k], ACTIVE, MASK1)) {
                    U[i ][j][k][_U] = UINIT;
                    U[i ][j][k][_V] = VINIT;
                    U[i ][j][k][_W] = WINIT;
                    P[i ][j][k]     = PINIT;
                    PP[i][j][k]     = 0;
                }
            }
        }
    }

    bc_u(U, UU, BU, X, KX, J, C, 0);
    bc_p(P, BP, U, UU, X, KX, J, C, 0);
}

void _clr_pp(void) {
    for (int i = I0; i <= I1; i ++) {
        for (int j = J0; j <= J1; j ++) {
            for (int k = K0; k <= K1; k ++) {
                if (f_see(F[i][j][k], ACTIVE, MASK1)) {
                    PP[i][j][k] = 0;
                }
            }
        }
    }
}

int main(void) {
    _init();

    _param_out();
    _var_out();
    _boundary_out();

    return 0;
}