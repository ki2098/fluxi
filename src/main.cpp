#include <stdio.h>
#include "tp.h"
#include "bc.h"
#include "var.h"
#include "ns.h"
#include "contra.h"
#include "solver.h"
#include "diver.h"
#include "turb.h"
#include "driver.h"

static void _param_out(void) {
    FILE *fo;
    fo = fopen("para.csv", "w+t");
    if (fo == NULL) {
        printf("\nERROR when opening file\n");
        fflush(stdout);
    }
    else {
        fprintf(fo, "x,y,z,active,fe,fn,ft,me,mn,mt,k1x1,k2x2,k3x3,j,g11,g22,g33\n");
        double x, y, z;
        int    flag, active;
        int    fe, fn, ft;
        int    me, mn, mt;
        double k1x1, k2x2, k3x3;
        double det;
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
                    fprintf(fo, "%.36lf,%.36lf,%.36lf,%d,%d,%d,%d,%d,%d,%d,%.36lf,%.36lf,%.36lf,%.36lf,%.36lf,%.36lf,%.36lf\n", x, y, z, active, fe, fn, ft, me, mn, mt, k1x1, k2x2, k3x3, det, g11, g22, g33);
                }
            }
        }
        fclose(fo);
    }
}

static void _var_out(char* fname) {
    #pragma acc data present(X, U, UT, P, SGS, DVR, DVA) 
    {
    // acc data starts
    #pragma acc update host(X, U, UT, P, SGS, DVR, DVA)
    
    FILE *fo;
    fo = fopen(fname, "w+t");
    if (fo == NULL) {
        printf("\nERROR when opening file\n");
        fflush(stdout);
    }
    else {
        fprintf(fo, "x,y,z,u,v,w,ut1,ut2,ut3,p,nue,dvr,dva\n");
        double x, y, z, u, v, w, p, nue, dvr, dva, ut1, ut2, ut3;
        for (int k = K0 - 1; k <= K1 + 1; k ++) {
            for (int j = J0 - 1; j <= J1 + 1; j ++) {
                for (int i = I0 - 1; i <= I1 + 1; i ++) {
                    x   =   X[i][j][k][_X];
                    y   =   X[i][j][k][_Y];
                    z   =   X[i][j][k][_Z];
                    u   =   U[i][j][k][_U];
                    v   =   U[i][j][k][_V];
                    w   =   U[i][j][k][_W];
                    ut1 =  UT[i][j][k][_E];
                    ut2 =  UT[i][j][k][_N];
                    ut3 =  UT[i][j][k][_T];
                    p   =   P[i][j][k];
                    nue = SGS[i][j][k];
                    dvr = DVR[i][j][k];
                    dva = DVA[i][j][k];
                    fprintf(fo, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", x, y, z, u, v, w, ut1, ut2, ut3, p, nue, dvr, dva);
                }
            }
        }
        fclose(fo);
    }

    // acc data ends
    }
}

static void _average_out(char* fname) {
    FILE *fo;
    fo = fopen(fname, "w+t");
    if (fo == NULL) {
        printf("\nERROR when opening file\n");
        fflush(stdout);
    }
    else {
        fprintf(fo, "x,y,z,u,v,w,p\n");
        double x, y, z, u, v, w, p;
        for (int k = K0; k <= K1; k ++) {
            for (int j = J0; j <= J1; j ++) {
                for (int i = I0; i <= I1; i ++) {
                    x =  X[i][j][k][_X];
                    y =  X[i][j][k][_Y];
                    z =  X[i][j][k][_Z];
                    u = UR[i][j][k][_U];
                    v = UR[i][j][k][_V];
                    w = UR[i][j][k][_W];
                    p = PR[i][j][k];
                    fprintf(fo, "%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", x, y, z, u, v, w, p);
                }
            }
        }
        fclose(fo);
    }
}

static void _boundary_out(void) {
    printf("nu du nv dv nw dw np dp\n");
    for (int i = 0; i <= BNUM; i ++) {
        unsigned int b  = B[i];
        unsigned int bu = f_see(b , _B_U, MASK2);
        unsigned int bv = f_see(b , _B_V, MASK2);
        unsigned int bw = f_see(b , _B_W, MASK2);
        unsigned int bp = f_see(b , _B_P, MASK2);
        unsigned int nu = f_see(bu, _NEUMANN, MASK1);
        unsigned int du = f_see(bu, _DIRICHLET, MASK1);
        unsigned int nv = f_see(bv, _NEUMANN, MASK1);
        unsigned int dv = f_see(bv, _DIRICHLET, MASK1);
        unsigned int nw = f_see(bw, _NEUMANN, MASK1);
        unsigned int dw = f_see(bw, _DIRICHLET, MASK1);
        unsigned int np = f_see(bp, _NEUMANN, MASK1);
        unsigned int dp = f_see(bp, _DIRICHLET, MASK1);
        printf("%u  %u  %u  %u  %u  %u  %u  %u\n", nu, du, nv, dv, nw, dw, np, dp);
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

static void _time_sum(void) {
    #pragma acc kernels loop independent collapse(3) present(U, UR, P, PR)
    for (int i = I0; i <= I1; i ++) {
        for (int j = J0 ; j <= J1; j ++) {
            for (int k = K0; k <= K1; k ++) {
                UR[i][j][k][_U] += U[i][j][k][_U];
                UR[i][j][k][_V] += U[i][j][k][_V];
                UR[i][j][k][_W] += U[i][j][k][_W];
                PR[i][j][k]     += P[i][j][k];
            }
        }
    }
}

static void _time_average(int steps) {
    for (int i = I0; i <= I1; i ++) {
        for (int j = J0 ; j <= J1; j ++) {
            for (int k = K0; k <= K1; k ++) {
                UR[i][j][k][_U] /= steps;
                UR[i][j][k][_V] /= steps;
                UR[i][j][k][_W] /= steps;
                PR[i][j][k]     /= steps;
            }
        }
    }
}

static void _init(void) {
    tp_x(X, KX, J, G);
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
    int    average_range  = 100000;
    int    step           = 0;
    double driver_p       = 0;
    double avg;
    double dva, dvr;
    double res;

    _init();
    _param_out();
    _boundary_out();

    printf(
        "%lf,%lf,%lf,%lf,%lf\n", 
        X[I0sq][J0sq][K1sq + 1][_Z] - Z1sq, 
        X0sq - X[I0sq - 1][J0sq][K0sq][_X], 
        X[I1sq + 1][J0sq][K0sq][_X] - X1sq,
        Y0sq - X[I0sq][J0sq - 1][K0sq][_Y],
        X[I0sq][J1sq + 1][K0sq][_Y] - Y1sq
    );

    #pragma acc enter data copyin(F, U, UA, UC, UU, UUA, UT, P, PD, DVR, DVA, SGS, UR, PR, X, KX, J, G, BU, BP)
    bc_u_periodic(U);
    bc_u_wall(F, U, UT, X);
    bc_p_driver(P, driver_p);
    bc_p_periodic(P);
    contra(F, U, UC, UU, BU, X, KX, J);
    turb_csm(F, U, BU, X, KX, J, SGS);

    sprintf(fname, "./data/var.csv.%d", n_file);
    _var_out(fname);
    n_file ++;

    for (step = 1; step <= NSTEP; step ++) {
        ns_pseudo_c(F, U, UA, UU, UT, BU, SGS, X, J, G);
        bc_u_periodic(UA);
        contra(F, UA, UC, UUA, BU, X, KX, J);
        diver(F, UUA, DVA, J, dva);

        driver_p_gradient(U, J, UINFLOW, driver_p);
        bc_p_driver(P, driver_p);

        iter_divergece = 0;
        do {
            iter_poisson = 0;
            do {
                solver_sor(F, P, BP, DVA, G, J, res);
                bc_p_driver(P, driver_p);
                bc_p_periodic(P);
                iter_poisson ++;
            } while (res > EPOI && iter_poisson < MAXIT);

            ns_correction_c(F, U, UA, P, BP, KX);
            bc_u_periodic(U);
            ns_correction_f(F, U, UU, UUA, BU, P, BP, X, KX, G);
            driver_monitor(U, J, avg);
            diver(F, UU, DVR, J, dvr);
            iter_divergece ++;
            printf("\rs(%6d,%6d):p(%4d,%13.10lf),d(%13.10lf,%13.10lf),u(%.6lf)", step, iter_divergece, iter_poisson, res, dva, dvr, avg);
            fflush(stdout);

            // if (iter_divergece >= 100) {
            //     goto END;
            // }
        } while (dvr > EDIV0);
        turb_csm(F, U, BU, X, KX, J, SGS);
        bc_u_wall(F, U, UT, X);

        _p_0_avg();
        bc_p_driver(P, driver_p);
        bc_p_periodic(P);

        if (step % int(1.0 / DT) == 0 || step == NSTEP) {
            sprintf(fname, "./data/var.csv.%d", n_file);
            _var_out(fname);
            n_file ++;
        }

        if (step > average_range) {
            _time_sum();
        }
    }
END:
    printf("\n");

    printf("NFile=%d\n", n_file);

    _var_out((char*)"./data/final.csv");

    #pragma acc exit data copyout(F, U, UA, UC, UU, UUA, UT, P, PD, DVR, DVA, SGS, UR, PR, X, KX, J, G, BU, BP)

    _time_average(step - average_range);
    _average_out((char*)"./data/time_average.csv");

    return 0;
}
