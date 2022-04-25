#include <stdio.h>
#include "tp.h"
#include "bc.h"
#include "var.h"
#include "ns.h"
#include "contra.h"
#include "solver.h"
#include "diver.h"

static void _param_out(void) {
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
        for (int k = K0 - 1; k <= K1 + 1; k ++) {
            for (int j = J0 - 1; j <= J1 + 1; j ++) {
                for (int i = I0 - 1; i <= I1 + 1; i ++) {
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
                    active = f_see(flag, _ACTIVE, MASK1);
                    fe     = f_see(flag, _F_E   , MASK8);
                    fn     = f_see(flag, _F_N   , MASK8);
                    ft     = f_see(flag, _F_T   , MASK8);
                    me     = f_see(flag, _M_E   , MASK1);
                    mn     = f_see(flag, _M_N   , MASK1);
                    mt     = f_see(flag, _M_T   , MASK1);
                    fprintf(fo, "%.36lf,%.36lf,%.36lf,%d,%d,%d,%d,%d,%d,%d,%.36lf,%.36lf,%.36lf,%.36lf,%.36lf,%.36lf,%.36lf,%.36lf,%.36lf,%.36lf,%.36lf,%.36lf,%.36lf\n", x, y, z, active, fe, fn, ft, me, mn, mt, k1x1, k2x2, k3x3, det, c1, c2, c3, c7, c8, c9, g11, g22, g33);
                }
            }
        }
        fclose(fo);
    }
}

static void _var_out(char* fname) {
    #pragma acc data present(X, U, UA, UP, UU, UUP, UUA, P, PP, DVR, DVP, DVA) 
    {
    // acc data starts
    #pragma acc update host(X, U, UA, UP, UU, UUP, UUA, P, PP, DVR, DVP, DVA)
    
    FILE *fo;
    fo = fopen(fname, "w+t");
    if (fo == NULL) {
        printf("\nERROR when opening file\n");
        fflush(stdout);
    }
    else {
        fprintf(fo, "x,y,z,u,v,w,ju,jv,jw,p,dvr,dva\n");
        double x, y, z, u, v, w, ju, jv, jw, p, dvr, dva;
        for (int k = K0; k <= K1; k ++) {
            for (int j = J0; j <= J1; j ++) {
                for (int i = I0; i <= I1; i ++) {
                    x   =   X[i][j][k][_X];
                    y   =   X[i][j][k][_Y];
                    z   =   X[i][j][k][_Z];
                    u   =   U[i][j][k][_U];
                    v   =   U[i][j][k][_V];
                    w   =   U[i][j][k][_W];
                    ju  =  UU[i][j][k][_U];
                    jv  =  UU[i][j][k][_V];
                    jw  =  UU[i][j][k][_W];
                    p   =   P[i][j][k];
                    dvr = DVR[i][j][k];
                    dva = DVA[i][j][k];
                    fprintf(fo, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", x, y, z, u, v, w, ju, jv, jw, p, dvr, dva);
                }
            }
        }
        fclose(fo);
    }

    // acc data ends
    }
}

static void _boundary_out(void) {
    printf("nu du nv dv nw dw np dp\n");
    for (int i = 0; i <= NB; i ++) {
        unsigned int b  = B[i];
        unsigned int bu = f_see(b, _BT_U, MASK2);
        unsigned int bv = f_see(b, _BT_V, MASK2);
        unsigned int bw = f_see(b, _BT_W, MASK2);
        unsigned int bp = f_see(b, _BT_P, MASK2);
        unsigned int nu = f_see(bu, _BT_N, MASK1);
        unsigned int du = f_see(bu, _BT_D, MASK1);
        unsigned int nv = f_see(bv, _BT_N, MASK1);
        unsigned int dv = f_see(bv, _BT_D, MASK1);
        unsigned int nw = f_see(bw, _BT_N, MASK1);
        unsigned int dw = f_see(bw, _BT_D, MASK1);
        unsigned int np = f_see(bp, _BT_N, MASK1);
        unsigned int dp = f_see(bp, _BT_D, MASK1);
        printf("%u  %u  %u  %u  %u  %u  %u  %u\n", nu, du, nv, dv, nw, dw, np, dp);
    }
}

static void _clr_pp(void) {
    #pragma acc kernels loop independent collapse(3) present(F, PP)
    for (int i = 0; i < NNX; i ++) {
        for (int j = 0; j < NNY; j ++) {
            for (int k = 0; k < NNZ; k ++) {
                if (f_see(F[i][j][k], _ACTIVE, MASK1)) {
                    PP[i][j][k] = 0;
                }
            }
        }
    }
}

static void _p_0_avg(void) {
    double sum = 0;
    int    cnt = 0;
    
    #pragma acc kernels loop independent collapse(3) reduction(+:sum, cnt) present(F, P) copy(sum, cnt)
    for (int i = I0; i <= I1; i ++) {
        for (int j = J0; j <= J1; j ++) {
            for (int k = K0; k <= K1; k ++) {
                if (f_see(F[i][j][k], _ACTIVE, MASK1)) {
                    sum += P[i][j][k];
                    cnt += 1;
                }
            }
        }
    }

    double avg = sum / cnt;

    #pragma acc kernels loop independent collapse(3) present(F, P) copyin(avg)
    for (int i = I0; i <= I1; i ++) {
        for (int j = J0; j <= J1; j ++) {
            for (int k = K0; k <= K1; k ++) {
                if (f_see(F[i][j][k], _ACTIVE, MASK1)) {
                    P[i][j][k] -= avg;
                }
            }
        }
    }
}

static void _init(void) {
    tp_x(X, KX, J, G, C);
    tp_f(F);
    
    for (int i = I0; i <= I1; i ++) {
        for (int j = J0; j <= J1; j ++) {
            for (int k = K0; k <= K1; k ++) {
                if (f_see(F[i][j][k], _ACTIVE, MASK1)) {
                    U[i ][j][k][_U] = UINIT;
                    U[i ][j][k][_V] = VINIT;
                    U[i ][j][k][_W] = WINIT;
                    P[i ][j][k]     = PINIT;
                }
            }
        }
    }

    bc_u_init(U, BU);
    bc_p_init(P, BP);
}

int main(void) {
    char   fname[128];
    int    n_file         = 0;
    int    iter_poisson   = 0;
    int    iter_divergece = 0;
    double dva, dvp, dvr;
    double res;

    _init();
    _param_out();
    _boundary_out();

    #pragma acc enter data copyin(F, U, UA, UP, UD, UC, UU, UUA, UUP, UUD, P, PD, PP, PPD, DVR, DVA, DVP, SGS, X, KX, J, G, C, BU, BP, BPP)
    contra(F, U, UC, UU, BU, X, KX, J);

    sprintf(fname, "./data/var.csv.%d", n_file);
    _var_out(fname);
    n_file ++;

    for (int step = 1; step <= NSTEP; step ++) {
        ns_pseudo_c(F, U, UA, UU, BU, SGS, KX, J, C);
        bc_u_outflow(U, UU, BU, X, KX);
        contra(F, UA, UC, UUA, BU, X, KX, J);
        diver(F, UUA, DVA, J, dva);

        // ns_correction_c(F, UP, UA, P, BP, KX);
        // ns_correction_f(F, UP, UUP, UUA, BU, P, BP, X, KX, G);
        // diver(F, UUP, DVP, J, dvp);

        // _clr_pp();
        iter_divergece = 0;
        do {
            iter_poisson = 0;
            do {
                // solver_sor(F, P, BP, DVP, C, res);
                solver_jacobi(F, P, PD, BP, DVA, C, res);
                bc_p_periodic(P);
                iter_poisson ++;
            } while (res > EPOI && iter_poisson < MAXIT);

            ns_correction_c(F, U, UA, P, BP, KX);
            ns_correction_f(F, U, UU, UUA, BU, P, BP, X, KX, G);
            diver(F, UU, DVR, J, dvr);
            iter_divergece ++;
            printf("\rs(%6d,%6d):p(%4d,%13.10lf),d(%13.10lf,%13.10lf)", step, iter_divergece, iter_poisson, res, dva, dvr);
            fflush(stdout);
        } while (dvr > EDIV);

        // ns_correction_p(P, PP);
        _p_0_avg();
        // bc_p_periodic(P);
        bc_u_periodic(U);

        if (step % int(0.5 / DT) == 0 || step == NSTEP /* || step <= 100 */) {
            sprintf(fname, "./data/var.csv.%d", n_file);
            _var_out(fname);
            n_file ++;
        }
    }
    printf("\n");

    printf("NFile=%d\n", n_file);

    #pragma acc exit data copyout(F, U, UA, UP, UD, UC, UU, UUA, UUP, UUD, P, PD, PP, PPD, DVR, DVA, DVP, SGS, X, KX, J, G, C, BU, BP, BPP)


    return 0;
}
