#include <math.h>
#include "turb.h"
#include "bc.h"

static void _pre_bc_eva(
    unsigned int  m,
    double        _0,
    double        _1,
    double       &ref,
    double       &dis
) {
    ref = (m)? _1 : _0;
    dis = 0.5 - m;
}

void turb_smagorinsky(
    unsigned int   F[NNX][NNY][NNZ],
    double         U[NNX][NNY][NNZ][3],
    double        BU[NNX][NNY][NNZ][3][3],
    double         X[NNX][NNY][NNZ][3],
    double        KX[NNX][NNY][NNZ][3],
    double         J[NNX][NNY][NNZ],
    double       SGS[NNX][NNY][NNZ]
) {
    #pragma acc kernels loop independent collapse(3) present(F, U, BU, X, KX, J, SGS) copyin(B)
    for (int i = I0; i <= I1; i ++) {
        for (int j = J0; j <= J1; j ++) {
            for (int k = K0; k <= K1; k ++) {
                if (f_see(F[i][j][k], _ACTIVE, MASK1)) {
                    unsigned int b12, b13;
                    unsigned int b22, b23;
                    unsigned int b32, b33;
                    unsigned int f12, f13;
                    unsigned int f22, f23;
                    unsigned int f32, f33;
                    unsigned int m12, m13;
                    unsigned int m22, m23;
                    unsigned int m32, m33;
                    double kx1, kx2, kx3;
                    double uc0;
                    double ue1, un1, ut1;
                    double uw1, us1, ub1;
                    double ref, dis;
                    double dux, duy, duz;
                    double dvx, dvy, dvz;
                    double dwx, dwy, dwz;
                    double det;
                    double a01, a02, a03;
                    double d01, d02, d03, d04, d05, d06;
                    double Del, D_u, l_0;

                    kx1 = KX[i    ][j    ][k    ][_X];
                    kx2 = KX[i    ][j    ][k    ][_Y];
                    kx3 = KX[i    ][j    ][k    ][_Z];
                    det =  J[i    ][j    ][k    ];
                    f12 =  F[i - 1][j    ][k    ];
                    f13 =  F[i    ][j    ][k    ];
                    f22 =  F[i    ][j - 1][k    ];
                    f23 =  F[i    ][j    ][k    ];
                    f32 =  F[i    ][j    ][k - 1];
                    f33 =  F[i    ][j    ][k    ];
                    m12 = f_see(f12, _M_E, MASK1);
                    m13 = f_see(f13, _M_E, MASK1);
                    m22 = f_see(f22, _M_N, MASK1);
                    m23 = f_see(f23, _M_N, MASK1);
                    m32 = f_see(f32, _M_T, MASK1);
                    m33 = f_see(f33, _M_T, MASK1);
                    f12 = f_see(f12, _F_E, MASK8);
                    f13 = f_see(f13, _F_E, MASK8);
                    f22 = f_see(f22, _F_N, MASK8);
                    f23 = f_see(f23, _F_N, MASK8);
                    f32 = f_see(f32, _F_T, MASK8);
                    f33 = f_see(f33, _F_T, MASK8);

                    uc0 = U[i    ][j    ][k    ][_U];
                    ue1 = U[i + 1][j    ][k    ][_U];
                    un1 = U[i    ][j + 1][k    ][_U];
                    ut1 = U[i    ][j    ][k + 1][_U];
                    uw1 = U[i - 1][j    ][k    ][_U];
                    us1 = U[i    ][j - 1][k    ][_U];
                    ub1 = U[i    ][j    ][k - 1][_U];
                    b12 = f_see(B[f12], _B_U, MASK2);
                    b13 = f_see(B[f13], _B_U, MASK2);
                    b22 = f_see(B[f22], _B_U, MASK2);
                    b23 = f_see(B[f23], _B_U, MASK2);
                    b32 = f_see(B[f32], _B_U, MASK2);
                    b33 = f_see(B[f33], _B_U, MASK2);
                    if (f13) {
                        _pre_bc_eva(m13, uc0, ue1, ref, dis);
                        ue1 = 2 * bc_eva(b13, ref, dis, BU[i][j][k][_E][_U]) - uc0;
                    }
                    if (f23) {
                        _pre_bc_eva(m23, uc0, un1, ref, dis);
                        un1 = 2 * bc_eva(b23, ref, dis, BU[i][j][k][_N][_U]) - uc0;
                    }
                    if (f33) {
                        _pre_bc_eva(m33, uc0, ut1, ref, dis);
                        ut1 = 2 * bc_eva(b33, ref, dis, BU[i][j][k][_T][_U]) - uc0;
                    }
                    if (f12) {
                        _pre_bc_eva(m12, uw1, uc0, ref, dis);
                        uw1 = 2 * bc_eva(b12, ref, dis, BU[i - 1][j][k][_E][_U]) - uc0;
                    }
                    if (f22) {
                        _pre_bc_eva(m22, us1, uc0, ref, dis);
                        us1 = 2 * bc_eva(b22, ref, dis, BU[i][j - 1][k][_N][_U]) - uc0;
                    }
                    if (f32) {
                        _pre_bc_eva(m32, ub1, uc0, ref, dis);
                        ub1 = 2 * bc_eva(b32, ref, dis, BU[i][j][k - 1][_T][_U]) - uc0;
                    }
                    dux = kx1 * 0.5 * (ue1 - uw1);
                    duy = kx2 * 0.5 * (un1 - us1);
                    duz = kx3 * 0.5 * (ut1 - ub1);

                    uc0 =   U[i    ][j    ][k    ][_V];
                    ue1 =   U[i + 1][j    ][k    ][_V];
                    un1 =   U[i    ][j + 1][k    ][_V];
                    ut1 =   U[i    ][j    ][k + 1][_V];
                    uw1 =   U[i - 1][j    ][k    ][_V];
                    us1 =   U[i    ][j - 1][k    ][_V];
                    ub1 =   U[i    ][j    ][k - 1][_V];
                    b12 = f_see(B[f12], _B_V, MASK2);
                    b13 = f_see(B[f13], _B_V, MASK2);
                    b22 = f_see(B[f22], _B_V, MASK2);
                    b23 = f_see(B[f23], _B_V, MASK2);
                    b32 = f_see(B[f32], _B_V, MASK2);
                    b33 = f_see(B[f33], _B_V, MASK2);
                    if (f13) {
                        _pre_bc_eva(m13, uc0, ue1, ref, dis);
                        ue1 = 2 * bc_eva(b13, ref, dis, BU[i][j][k][_E][_V]) - uc0;
                    }
                    if (f23) {
                        _pre_bc_eva(m23, uc0, un1, ref, dis);
                        un1 = 2 * bc_eva(b23, ref, dis, BU[i][j][k][_N][_V]) - uc0;
                    }
                    if (f33) {
                        _pre_bc_eva(m33, uc0, ut1, ref, dis);
                        ut1 = 2 * bc_eva(b33, ref, dis, BU[i][j][k][_T][_V]) - uc0;
                    }
                    if (f12) {
                        _pre_bc_eva(m12, uw1, uc0, ref, dis);
                        uw1 = 2 * bc_eva(b12, ref, dis, BU[i - 1][j][k][_E][_V]) - uc0;
                    }
                    if (f22) {
                        _pre_bc_eva(m22, us1, uc0, ref, dis);
                        us1 = 2 * bc_eva(b22, ref, dis, BU[i][j - 1][k][_N][_V]) - uc0;
                    }
                    if (f32) {
                        _pre_bc_eva(m32, ub1, uc0, ref, dis);
                        ub1 = 2 * bc_eva(b32, ref, dis, BU[i][j][k - 1][_T][_V]) - uc0;
                    }
                    dvx = kx1 * 0.5 * (ue1 - uw1);
                    dvy = kx2 * 0.5 * (un1 - us1);
                    dvz = kx3 * 0.5 * (ut1 - ub1);

                    uc0 =   U[i    ][j    ][k    ][_W];
                    ue1 =   U[i + 1][j    ][k    ][_W];
                    un1 =   U[i    ][j + 1][k    ][_W];
                    ut1 =   U[i    ][j    ][k + 1][_W];
                    uw1 =   U[i - 1][j    ][k    ][_W];
                    us1 =   U[i    ][j - 1][k    ][_W];
                    ub1 =   U[i    ][j    ][k - 1][_W];
                    b12 = f_see(B[f12], _B_W, MASK2);
                    b13 = f_see(B[f13], _B_W, MASK2);
                    b22 = f_see(B[f22], _B_W, MASK2);
                    b23 = f_see(B[f23], _B_W, MASK2);
                    b32 = f_see(B[f32], _B_W, MASK2);
                    b33 = f_see(B[f33], _B_W, MASK2);
                    if (f13) {
                        _pre_bc_eva(m13, uc0, ue1, ref, dis);
                        ue1 = 2 * bc_eva(b13, ref, dis, BU[i][j][k][_E][_W]) - uc0;
                    }
                    if (f23) {
                        _pre_bc_eva(m23, uc0, un1, ref, dis);
                        un1 = 2 * bc_eva(b23, ref, dis, BU[i][j][k][_N][_W]) - uc0;
                    }
                    if (f33) {
                        _pre_bc_eva(m33, uc0, ut1, ref, dis);
                        ut1 = 2 * bc_eva(b33, ref, dis, BU[i][j][k][_T][_W]) - uc0;
                    }
                    if (f12) {
                        _pre_bc_eva(m12, uw1, uc0, ref, dis);
                        uw1 = 2 * bc_eva(b12, ref, dis, BU[i - 1][j][k][_E][_W]) - uc0;
                    }
                    if (f22) {
                        _pre_bc_eva(m22, us1, uc0, ref, dis);
                        us1 = 2 * bc_eva(b22, ref, dis, BU[i][j - 1][k][_N][_W]) - uc0;
                    }
                    if (f32) {
                        _pre_bc_eva(m32, ub1, uc0, ref, dis);
                        ub1 = 2 * bc_eva(b32, ref, dis, BU[i][j][k - 1][_T][_W]) - uc0;
                    }
                    dwx = kx1 * 0.5 * (ue1 - uw1);
                    dwy = kx2 * 0.5 * (un1 - us1);
                    dwz = kx3 * 0.5 * (ut1 - ub1);

                    
                    a01 = U[i][j][K1][_U] / (X[i][j][K1][_Z] - Z0);
                    a02 = (X[i][j][k][_Z] - Z0) * sqrt(RE * fabs(a01));
                    a03 = 1.0 - exp(- a02 / 25.0);
                    d01 = 2 * dux * dux;
                    d02 = 2 * dvy * dvy;
                    d03 = 2 * dwz * dwz;
                    d04 = (dwy + dvz) * (dwy + dvz);
                    d05 = (duz + dwx) * (duz + dwx);
                    d06 = (duy + dvx) * (duy + dvx);
                    D_u = sqrt(d01 + d02 + d03 + d04 + d05 + d06);
                    Del = cbrt(det);
                    l_0 = C_s * a03 * Del;

                    SGS[i][j][k] = l_0 * l_0 * D_u;
                }
            }
        }
    }
}
