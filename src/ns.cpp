#include <math.h>
#include <stdio.h>
#include "ns.h"
#include "util.h"
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

static double _ffvc_muscl(
    double uc0,
    double ue1,
    double ue2,
    double un1,
    double un2,
    double ut1,
    double ut2,
    double uw1,
    double uw2,
    double us1,
    double us2,
    double ub1,
    double ub2,
    double ufe,
    double vfn,
    double wft,
    double ufw,
    double vfs,
    double wfb,
    double det,
    int    epe,
    int    epn,
    int    ept,
    int    epw,
    int    eps,
    int    epb
) {
    double d1, d2, d3, d4;
    double s1, s2, s3, s4;
    double g1, g2, g3, g4, g5, g6;
    double u11, u10, u01, u00;
    double f1, f0;
    double adv;

    double kappa = 1.0 / 3.0;
    double beta  = (3.0 - kappa) / (1.0 - kappa);
    double m1    = 1.0 - kappa;
    double m2    = 1.0 + kappa;

    d4 = ue2 - ue1;
    d3 = ue1 - uc0;
    d2 = uc0 - uw1;
    d1 = uw1 - uw2;
    s4 = copysign(1.0, d4);
    s3 = copysign(1.0, d3);
    s2 = copysign(1.0, d2);
    s1 = copysign(1.0, d1);
    g6  = s4 * max(0.0, min(fabs(d4), s4 * beta * d3));
    g5  = s3 * max(0.0, min(fabs(d3), s3 * beta * d4));
    g4  = s3 * max(0.0, min(fabs(d3), s3 * beta * d2));
    g3  = s2 * max(0.0, min(fabs(d2), s2 * beta * d3));
    g2  = s2 * max(0.0, min(fabs(d2), s2 * beta * d1));
    g1  = s1 * max(0.0, min(fabs(d1), s1 * beta * d2));
    u11 = ue1 - 0.25 * (m1 * g6 + m2 * g5) * epe;
    u10 = uc0 + 0.25 * (m1 * g3 + m2 * g4);
    u01 = uc0 - 0.25 * (m1 * g4 + m2 * g3);
    u00 = uw1 + 0.25 * (m1 * g1 + m2 * g2) * epw;
    f1  = 0.5 * (ufe * (u11 + u10) - fabs(ufe) * (u11 - u10));
    f0  = 0.5 * (ufw * (u01 + u00) - fabs(ufw) * (u01 - u00));
    adv = (f1 - f0) / det;

    d4  = un2 - un1;
    d3  = un1 - uc0;
    d2  = uc0 - us1;
    d1  = us1 - us2;
    s4  = copysign(1.0, d4);
    s3  = copysign(1.0, d3);
    s2  = copysign(1.0, d2);
    s1  = copysign(1.0, d1);
    g6  = s4 * max(0.0, min(fabs(d4), s4 * beta * d3));
    g5  = s3 * max(0.0, min(fabs(d3), s3 * beta * d4));
    g4  = s3 * max(0.0, min(fabs(d3), s3 * beta * d2));
    g3  = s2 * max(0.0, min(fabs(d2), s2 * beta * d3));
    g2  = s2 * max(0.0, min(fabs(d2), s2 * beta * d1));
    g1  = s1 * max(0.0, min(fabs(d1), s1 * beta * d2));
    u11 = un1 - 0.25 * (m1 * g6 + m2 * g5) * epn;
    u10 = uc0 + 0.25 * (m1 * g3 + m2 * g4);
    u01 = uc0 - 0.25 * (m1 * g4 + m2 * g3);
    u00 = us1 + 0.25 * (m1 * g1 + m2 * g2) * eps;
    f1  = 0.5 * (vfn * (u11 + u10) - fabs(vfn) * (u11 - u10));
    f0  = 0.5 * (vfs * (u01 + u00) - fabs(vfs) * (u01 - u00));
    adv += (f1 - f0) / det;

    d4  = ut2 - ut1;
    d3  = ut1 - uc0;
    d2  = uc0 - ub1;
    d1  = ub1 - ub2;
    s4  = copysign(1.0, d4);
    s3  = copysign(1.0, d3);
    s2  = copysign(1.0, d2);
    s1  = copysign(1.0, d1);
    g6  = s4 * max(0.0, min(fabs(d4), s4 * beta * d3));
    g5  = s3 * max(0.0, min(fabs(d3), s3 * beta * d4));
    g4  = s3 * max(0.0, min(fabs(d3), s3 * beta * d2));
    g3  = s2 * max(0.0, min(fabs(d2), s2 * beta * d3));
    g2  = s2 * max(0.0, min(fabs(d2), s2 * beta * d1));
    g1  = s1 * max(0.0, min(fabs(d1), s1 * beta * d2));
    u11 = ut1 - 0.25 * (m1 * g6 + m2 * g5) * ept;
    u10 = uc0 + 0.25 * (m1 * g3 + m2 * g4);
    u01 = uc0 - 0.25 * (m1 * g4 + m2 * g3);
    u00 = ub1 + 0.25 * (m1 * g1 + m2 * g2) * epb;
    f1  = 0.5 * (wft * (u11 + u10) - fabs(wft) * (u11 - u10));
    f0  = 0.5 * (wfb * (u01 + u00) - fabs(wfb) * (u01 - u00));
    adv += (f1 - f0) / det;

    return adv;
}

static double _vis(
    unsigned int f13,
    unsigned int f23,
    unsigned int f33,
    unsigned int f12,
    unsigned int f22,
    unsigned int f32,
    double       mag,
    double       uc0,
    double       ue1,
    double       un1,
    double       ut1,
    double       uw1,
    double       us1,
    double       ub1,
    double       ute,
    double       utn,
    double       utt,
    double       utw,
    double       uts,
    double       utb,
    double       de1,
    double       dn1,
    double       dt1,
    double       dw1,
    double       ds1,
    double       db1,
    double       nc0,
    double       ne1,
    double       nn1,
    double       nt1,
    double       nw1,
    double       ns1,
    double       nb1,
    double       det,
    double       g1c,
    double       g2c,
    double       g3c,
    double       g1e,
    double       g2n,
    double       g3t,
    double       g1w,
    double       g2s,
    double       g3b
) {
    double ffe, ffn, fft, ffw, ffs, ffb;
    ffe = (RI + 0.5 * (nc0 + ne1)) * 0.5 * (g1e + g1c) * (ue1 - uc0);
    ffn = (RI + 0.5 * (nc0 + nn1)) * 0.5 * (g2n + g2c) * (un1 - uc0);
    fft = (RI + 0.5 * (nc0 + nt1)) * 0.5 * (g3t + g3c) * (ut1 - uc0);
    ffw = (RI + 0.5 * (nc0 + nw1)) * 0.5 * (g1c + g1w) * (uc0 - uw1);
    ffs = (RI + 0.5 * (nc0 + ns1)) * 0.5 * (g2c + g2s) * (uc0 - us1);
    ffb = (RI + 0.5 * (nc0 + nb1)) * 0.5 * (g3c + g3b) * (uc0 - ub1);
    if (f13 == WALL) {
        if (mag < 1e-6) {
            ffe = 0;
        } else {
            double uti = ute * uc0 / mag;
            ffe = copysign(1.0, uc0) * 0.5 * (g1e + g1c) * de1 * uti * uti;
        }
    }
    if (f23 == WALL) {
        if (mag < 1e-6) {
            ffn = 0;
        } else {
            double uti = utn * uc0 / mag;
            ffn = copysign(1.0, uc0) * 0.5 * (g2n + g2c) * dn1 * uti * uti;
        }
    }
    if (f33 == WALL) {
        if (mag < 1e-6) {
            fft = 0;
        } else {
            double uti = utt * uc0 / mag;
            fft = copysign(1.0, uc0) * 0.5 * (g3t + g3c) * dt1 * uti * uti;
        }
    }
    if (f12 == WALL) {
        if (mag < 1e-6) {
            ffw = 0;
        } else {
            double uti = utw * uc0 / mag;
            ffw = copysign(1.0, uc0) * 0.5 * (g1c + g1w) * dw1 * uti * uti;
        }
    }
    if (f22 == WALL) {
        if (mag < 1e-6) {
            ffs = 0;
        } else {
            double uti = uts * uc0 / mag;
            ffs = copysign(1.0, uc0) * 0.5 * (g2c + g2s) * ds1 * uti * uti;
        }
    }
    if (f32 == WALL) {
        if (mag < 1e-6) {
            ffb = 0;
        } else {
            double uti = utb * uc0 / mag;
            ffb = copysign(1.0, uc0) * 0.5 * (g3c + g3b) * db1 * uti * uti;
        }
    }
    return (ffe + ffn + fft - ffw - ffs - ffb) / det;
}

void ns_pseudo_c(
    unsigned int   F[NNX][NNY][NNZ],
    double         U[NNX][NNY][NNZ][3],
    double        UA[NNX][NNY][NNZ][3],
    double        UU[NNX][NNY][NNZ][3],
    double        UT[NNX][NNY][NNZ][3],
    double        BU[NNX][NNY][NNZ][3][3],
    double       SGS[NNX][NNY][NNZ],
    double         X[NNX][NNY][NNZ][3],
    double         J[NNX][NNY][NNZ],
    double         G[NNX][NNY][NNZ][3]
) {
    #pragma acc kernels loop independent collapse(3) present(F, U, UA, UU, UT, BU, SGS, X, J, G) copyin(B)
    for (int i = I0; i <= I1; i ++) {
        for (int j = J0; j <= J1; j ++) {
            for (int k = K0; k <= K1; k ++) {
                if (f_see(F[i][j][k], _ACTIVE, MASK1)) {
                    double uc0;
                    double ue1, un1, ut1;
                    double ue2, un2, ut2;
                    double uw1, us1, ub1;
                    double uw2, us2, ub2;
                    double ufe, vfn, wft;
                    double ufw, vfs, wfb;
                    double ute, utn, utt;
                    double utw, uts, utb;
                    double de1, dn1, dt1;
                    double dw1, ds1, db1;
                    double nc0;
                    double ne1, nn1, nt1;
                    double nw1, ns1, nb1;
                    double det;
                    double ref, dis;
                    double ad1, ad2, ad3;
                    double vi1, vi2, vi3;
                    double xc0, yc0, zc0;
                    double xe1, yn1, zt1;
                    double xw1, ys1, zb1;
                    double g1c, g2c, g3c;
                    double g1e, g2n, g3t;
                    double g1w, g2s, g3b;
                    unsigned int f11, f12, f13, f14;
                    unsigned int f21, f22, f23, f24;
                    unsigned int f31, f32, f33, f34;
                    unsigned int b11, b12, b13, b14;
                    unsigned int b21, b22, b23, b24;
                    unsigned int b31, b32, b33, b34;
                    unsigned int m11, m12, m13, m14;
                    unsigned int m21, m22, m23, m24;
                    unsigned int m31, m32, m33, m34;
                    int epe, epn, ept;
                    int epw, eps, epb;
                    double u1, u2, u3, mag;
                    u1  = U[i][j][k][_U];
                    u2  = U[i][j][k][_V];
                    u3  = U[i][j][k][_W];
                    mag = sqrt(u1 * u1 + u2 * u2 + u3 * u3);

                    ufe =  UU[i    ][j    ][k    ][_U];
                    vfn =  UU[i    ][j    ][k    ][_V];
                    wft =  UU[i    ][j    ][k    ][_W];
                    ufw =  UU[i - 1][j    ][k    ][_U];
                    vfs =  UU[i    ][j - 1][k    ][_V];
                    wfb =  UU[i    ][j    ][k - 1][_W];
                    ute =  UT[i    ][j    ][k    ][_E];
                    utn =  UT[i    ][j    ][k    ][_N];
                    utt =  UT[i    ][j    ][k    ][_T];
                    utw =  UT[i - 1][j    ][k    ][_E];
                    uts =  UT[i    ][j - 1][k    ][_N];
                    utb =  UT[i    ][j    ][k - 1][_T];
                    det =   J[i    ][j    ][k    ];

                    g1c =   G[i    ][j    ][k    ][_X];
                    g2c =   G[i    ][j    ][k    ][_Y];
                    g3c =   G[i    ][j    ][k    ][_Z];
                    g1e =   G[i + 1][j    ][k    ][_X];
                    g2n =   G[i    ][j + 1][k    ][_Y];
                    g3t =   G[i    ][j    ][k + 1][_Z];
                    g1w =   G[i - 1][j    ][k    ][_X];
                    g2s =   G[i    ][j - 1][k    ][_Y];
                    g3b =   G[i    ][j    ][k - 1][_Z];
                    xc0 =   X[i    ][j    ][k    ][_X];
                    yc0 =   X[i    ][j    ][k    ][_Y];
                    zc0 =   X[i    ][j    ][k    ][_Z];
                    xe1 =   X[i + 1][j    ][k    ][_X];
                    yn1 =   X[i    ][j + 1][k    ][_Y];
                    zt1 =   X[i    ][j    ][k + 1][_Z];
                    xw1 =   X[i - 1][j    ][k    ][_X];
                    ys1 =   X[i    ][j - 1][k    ][_Y];
                    zb1 =   X[i    ][j    ][k - 1][_Z];
                    f11 =   F[i - 2][j    ][k    ];
                    f12 =   F[i - 1][j    ][k    ];
                    f13 =   F[i    ][j    ][k    ];
                    f14 =   F[i + 1][j    ][k    ];
                    f21 =   F[i    ][j - 2][k    ];
                    f22 =   F[i    ][j - 1][k    ];
                    f23 =   F[i    ][j    ][k    ];
                    f24 =   F[i    ][j + 1][k    ];
                    f31 =   F[i    ][j    ][k - 2];
                    f32 =   F[i    ][j    ][k - 1];
                    f33 =   F[i    ][j    ][k    ];
                    f34 =   F[i    ][j    ][k + 1];
                    m11 = f_see(f11, _M_E, MASK1);
                    m12 = f_see(f12, _M_E, MASK1);
                    m13 = f_see(f13, _M_E, MASK1);
                    m14 = f_see(f14, _M_E, MASK1);
                    m21 = f_see(f21, _M_N, MASK1);
                    m22 = f_see(f22, _M_N, MASK1);
                    m23 = f_see(f23, _M_N, MASK1);
                    m24 = f_see(f24, _M_N, MASK1);
                    m31 = f_see(f31, _M_T, MASK1);
                    m32 = f_see(f32, _M_T, MASK1);
                    m33 = f_see(f33, _M_T, MASK1);
                    m34 = f_see(f34, _M_T, MASK1);
                    f11 = f_see(f11, _F_E, MASK8);
                    f12 = f_see(f12, _F_E, MASK8);
                    f13 = f_see(f13, _F_E, MASK8);
                    f14 = f_see(f14, _F_E, MASK8);
                    f21 = f_see(f21, _F_N, MASK8);
                    f22 = f_see(f22, _F_N, MASK8);
                    f23 = f_see(f23, _F_N, MASK8);
                    f24 = f_see(f24, _F_N, MASK8);
                    f31 = f_see(f31, _F_T, MASK8);
                    f32 = f_see(f32, _F_T, MASK8);
                    f33 = f_see(f33, _F_T, MASK8);
                    f34 = f_see(f34, _F_T, MASK8);
                    epe = (f13)? 0 : 1;
                    epn = (f23)? 0 : 1;
                    ept = (f33)? 0 : 1;
                    epw = (f12)? 0 : 1;
                    eps = (f22)? 0 : 1;
                    epb = (f32)? 0 : 1;
                    de1 = xc0 - xe1;
                    dn1 = yc0 - yn1;
                    dt1 = zc0 - zt1;
                    dw1 = xc0 - xw1;
                    ds1 = yc0 - ys1;
                    db1 = zc0 - zb1;

                    nc0 = SGS[i    ][j    ][k    ];
                    ne1 = SGS[i + 1][j    ][k    ];
                    nn1 = SGS[i    ][j + 1][k    ];
                    nt1 = SGS[i    ][j    ][k + 1];
                    nw1 = SGS[i - 1][j    ][k    ];
                    ns1 = SGS[i    ][j - 1][k    ];
                    nb1 = SGS[i    ][j    ][k - 1];
                    b12 = f_see(B[f12], _B_Nt, MASK2);
                    b13 = f_see(B[f13], _B_Nt, MASK2);
                    b22 = f_see(B[f22], _B_Nt, MASK2);
                    b23 = f_see(B[f23], _B_Nt, MASK2);
                    b32 = f_see(B[f32], _B_Nt, MASK2);
                    b33 = f_see(B[f33], _B_Nt, MASK2);
                    if (f13) {
                        _pre_bc_eva(m13, nc0, ne1, ref, dis);
                        ne1 = 2 * bc_eva(b13, ref, dis, Nt_BOUNDARY_VALUE) - nc0;
                    }
                    if (f23) {
                        _pre_bc_eva(m23, nc0, nn1, ref, dis);
                        nn1 = 2 * bc_eva(b23, ref, dis, Nt_BOUNDARY_VALUE) - nc0;
                    }
                    if (f33) {
                        _pre_bc_eva(m33, nc0, nt1, ref, dis);
                        nt1 = 2 * bc_eva(b33, ref, dis, Nt_BOUNDARY_VALUE) - nc0;
                    }
                    if (f12) {
                        _pre_bc_eva(m12, nw1, nc0, ref, dis);
                        nw1 = 2 * bc_eva(b12, ref, dis, Nt_BOUNDARY_VALUE) - nc0;
                    }
                    if (f22) {
                        _pre_bc_eva(m22, ns1, nc0, ref, dis);
                        ns1 = 2 * bc_eva(b22, ref, dis, Nt_BOUNDARY_VALUE) - nc0;
                    }
                    if (f32) {
                        _pre_bc_eva(m32, nb1, nc0, ref, dis);
                        nb1 = 2 * bc_eva(b32, ref, dis, Nt_BOUNDARY_VALUE) - nc0;
                    }

                    uc0 = u1;
                    ue1 = U[i + 1][j    ][k    ][_U];
                    ue2 = U[i + 2][j    ][k    ][_U];
                    un1 = U[i    ][j + 1][k    ][_U];
                    un2 = U[i    ][j + 2][k    ][_U];
                    ut1 = U[i    ][j    ][k + 1][_U];
                    ut2 = U[i    ][j    ][k + 2][_U];
                    uw1 = U[i - 1][j    ][k    ][_U];
                    uw2 = U[i - 2][j    ][k    ][_U];
                    us1 = U[i    ][j - 1][k    ][_U];
                    us2 = U[i    ][j - 2][k    ][_U];
                    ub1 = U[i    ][j    ][k - 1][_U];
                    ub2 = U[i    ][j    ][k - 2][_U];
                    b11 = f_see(B[f11], _B_U, MASK2);
                    b12 = f_see(B[f12], _B_U, MASK2);
                    b13 = f_see(B[f13], _B_U, MASK2);
                    b14 = f_see(B[f14], _B_U, MASK2);
                    b21 = f_see(B[f21], _B_U, MASK2);
                    b22 = f_see(B[f22], _B_U, MASK2);
                    b23 = f_see(B[f23], _B_U, MASK2);
                    b24 = f_see(B[f24], _B_U, MASK2);
                    b31 = f_see(B[f31], _B_U, MASK2);
                    b32 = f_see(B[f32], _B_U, MASK2);
                    b33 = f_see(B[f33], _B_U, MASK2);
                    b34 = f_see(B[f34], _B_U, MASK2);
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
                    if (f14) {
                        _pre_bc_eva(m14, ue1, ue2, ref, dis);
                        ue2 = 2 * bc_eva(b14, ref, dis, BU[i + 1][j][k][_E][_U]) - ue1;
                    }
                    if (f24) {
                        _pre_bc_eva(m24, un1, un2, ref, dis);
                        un2 = 2 * bc_eva(b24, ref, dis, BU[i][j + 1][k][_N][_U]) - un1;
                    }
                    if (f34) {
                        _pre_bc_eva(m34, ut1, ut2, ref, dis);
                        ut2 = 2 * bc_eva(b34, ref, dis, BU[i][j][k + 1][_T][_U]) - ut1;
                    }
                    if (f11) {
                        _pre_bc_eva(m11, uw2, uw1, ref, dis);
                        uw2 = 2 * bc_eva(b11, ref, dis, BU[i - 2][j][k][_E][_U]) - uw1;
                    }
                    if (f21) {
                        _pre_bc_eva(m21, us2, us1, ref, dis);
                        us2 = 2 * bc_eva(b21, ref, dis, BU[i][j - 2][k][_N][_U]) - us1;
                    }
                    if (f31) {
                        _pre_bc_eva(m31, ub2, ub1, ref, dis);
                        ub2 = 2 * bc_eva(b31, ref, dis, BU[i][j][k - 2][_T][_U]) - ub1;
                    }
                    ad1 = _ffvc_muscl(uc0, ue1, ue2, un1, un2, ut1, ut2, uw1, uw2, us1, us2, ub1, ub2, ufe, vfn, wft, ufw, vfs, wfb, det, epe, epn, ept, epw, eps, epb);
                    vi1 = _vis(f13, f23, f33, f12, f22, f32, mag, uc0, ue1, un1, ut1, uw1, us1, ub1, ute, utn, utt, utw, uts, utb, de1, dn1, dt1, dw1, ds1, db1, nc0, ne1, nn1, nt1, nw1, ns1, nb1, det, g1c, g2c, g3c, g1e, g2n, g3t, g1w, g2s, g3b);

                    uc0 =   u2;
                    ue1 =   U[i + 1][j    ][k    ][_V];
                    ue2 =   U[i + 2][j    ][k    ][_V];
                    un1 =   U[i    ][j + 1][k    ][_V];
                    un2 =   U[i    ][j + 2][k    ][_V];
                    ut1 =   U[i    ][j    ][k + 1][_V];
                    ut2 =   U[i    ][j    ][k + 2][_V];
                    uw1 =   U[i - 1][j    ][k    ][_V];
                    uw2 =   U[i - 2][j    ][k    ][_V];
                    us1 =   U[i    ][j - 1][k    ][_V];
                    us2 =   U[i    ][j - 2][k    ][_V];
                    ub1 =   U[i    ][j    ][k - 1][_V];
                    ub2 =   U[i    ][j    ][k - 2][_V];
                    b11 = f_see(B[f11], _B_V, MASK2);
                    b12 = f_see(B[f12], _B_V, MASK2);
                    b13 = f_see(B[f13], _B_V, MASK2);
                    b14 = f_see(B[f14], _B_V, MASK2);
                    b21 = f_see(B[f21], _B_V, MASK2);
                    b22 = f_see(B[f22], _B_V, MASK2);
                    b23 = f_see(B[f23], _B_V, MASK2);
                    b24 = f_see(B[f24], _B_V, MASK2);
                    b31 = f_see(B[f31], _B_V, MASK2);
                    b32 = f_see(B[f32], _B_V, MASK2);
                    b33 = f_see(B[f33], _B_V, MASK2);
                    b34 = f_see(B[f34], _B_V, MASK2);
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
                    if (f14) {
                        _pre_bc_eva(m14, ue1, ue2, ref, dis);
                        ue2 = 2 * bc_eva(b14, ref, dis, BU[i + 1][j][k][_E][_V]) - ue1;
                    }
                    if (f24) {
                        _pre_bc_eva(m24, un1, un2, ref, dis);
                        un2 = 2 * bc_eva(b24, ref, dis, BU[i][j + 1][k][_N][_V]) - un1;
                    }
                    if (f34) {
                        _pre_bc_eva(m34, ut1, ut2, ref, dis);
                        ut2 = 2 * bc_eva(b34, ref, dis, BU[i][j][k + 1][_T][_V]) - ut1;
                    }
                    if (f11) {
                        _pre_bc_eva(m11, uw2, uw1, ref, dis);
                        uw2 = 2 * bc_eva(b11, ref, dis, BU[i - 2][j][k][_E][_V]) - uw1;
                    }
                    if (f21) {
                        _pre_bc_eva(m21, us2, us1, ref, dis);
                        us2 = 2 * bc_eva(b21, ref, dis, BU[i][j - 2][k][_N][_V]) - us1;
                    }
                    if (f31) {
                        _pre_bc_eva(m31, ub2, ub1, ref, dis);
                        ub2 = 2 * bc_eva(b31, ref, dis, BU[i][j][k - 2][_T][_V]) - ub1;
                    }
                    ad2 = _ffvc_muscl(uc0, ue1, ue2, un1, un2, ut1, ut2, uw1, uw2, us1, us2, ub1, ub2, ufe, vfn, wft, ufw, vfs, wfb, det, epe, epn, ept, epw, eps, epb);
                    vi2 = _vis(f13, f23, f33, f12, f22, f32, mag, uc0, ue1, un1, ut1, uw1, us1, ub1, ute, utn, utt, utw, uts, utb, de1, dn1, dt1, dw1, ds1, db1, nc0, ne1, nn1, nt1, nw1, ns1, nb1, det, g1c, g2c, g3c, g1e, g2n, g3t, g1w, g2s, g3b);

                    uc0 =   u3;
                    ue1 =   U[i + 1][j    ][k    ][_W];
                    ue2 =   U[i + 2][j    ][k    ][_W];
                    un1 =   U[i    ][j + 1][k    ][_W];
                    un2 =   U[i    ][j + 2][k    ][_W];
                    ut1 =   U[i    ][j    ][k + 1][_W];
                    ut2 =   U[i    ][j    ][k + 2][_W];
                    uw1 =   U[i - 1][j    ][k    ][_W];
                    uw2 =   U[i - 2][j    ][k    ][_W];
                    us1 =   U[i    ][j - 1][k    ][_W];
                    us2 =   U[i    ][j - 2][k    ][_W];
                    ub1 =   U[i    ][j    ][k - 1][_W];
                    ub2 =   U[i    ][j    ][k - 2][_W];
                    b11 = f_see(B[f11], _B_W, MASK2);
                    b12 = f_see(B[f12], _B_W, MASK2);
                    b13 = f_see(B[f13], _B_W, MASK2);
                    b14 = f_see(B[f14], _B_W, MASK2);
                    b21 = f_see(B[f21], _B_W, MASK2);
                    b22 = f_see(B[f22], _B_W, MASK2);
                    b23 = f_see(B[f23], _B_W, MASK2);
                    b24 = f_see(B[f24], _B_W, MASK2);
                    b31 = f_see(B[f31], _B_W, MASK2);
                    b32 = f_see(B[f32], _B_W, MASK2);
                    b33 = f_see(B[f33], _B_W, MASK2);
                    b34 = f_see(B[f34], _B_W, MASK2);
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
                    if (f14) {
                        _pre_bc_eva(m14, ue1, ue2, ref, dis);
                        ue2 = 2 * bc_eva(b14, ref, dis, BU[i + 1][j][k][_E][_W]) - ue1;
                    }
                    if (f24) {
                        _pre_bc_eva(m24, un1, un2, ref, dis);
                        un2 = 2 * bc_eva(b24, ref, dis, BU[i][j + 1][k][_N][_W]) - un1;
                    }
                    if (f34) {
                        _pre_bc_eva(m34, ut1, ut2, ref, dis);
                        ut2 = 2 * bc_eva(b34, ref, dis, BU[i][j][k + 1][_T][_W]) - ut1;
                    }
                    if (f11) {
                        _pre_bc_eva(m11, uw2, uw1, ref, dis);
                        uw2 = 2 * bc_eva(b11, ref, dis, BU[i - 2][j][k][_E][_W]) - uw1;
                    }
                    if (f21) {
                        _pre_bc_eva(m21, us2, us1, ref, dis);
                        us2 = 2 * bc_eva(b21, ref, dis, BU[i][j - 2][k][_N][_W]) - us1;
                    }
                    if (f31) {
                        _pre_bc_eva(m31, ub2, ub1, ref, dis);
                        ub2 = 2 * bc_eva(b31, ref, dis, BU[i][j][k - 2][_T][_W]) - ub1;
                    }
                    ad3 = _ffvc_muscl(uc0, ue1, ue2, un1, un2, ut1, ut2, uw1, uw2, us1, us2, ub1, ub2, ufe, vfn, wft, ufw, vfs, wfb, det, epe, epn, ept, epw, eps, epb);
                    vi3 = _vis(f13, f23, f33, f12, f22, f32, mag, uc0, ue1, un1, ut1, uw1, us1, ub1, ute, utn, utt, utw, uts, utb, de1, dn1, dt1, dw1, ds1, db1, nc0, ne1, nn1, nt1, nw1, ns1, nb1, det, g1c, g2c, g3c, g1e, g2n, g3t, g1w, g2s, g3b);

                    UA[i][j][k][_U] = U[i][j][k][_U] + DT * (- ad1 + vi1);
                    UA[i][j][k][_V] = U[i][j][k][_V] + DT * (- ad2 + vi2);
                    UA[i][j][k][_W] = U[i][j][k][_W] + DT * (- ad3 + vi3);
                }
            }
        }
    }
}

void ns_correction_c(
    unsigned int   F[NNX][NNY][NNZ],
    double         U[NNX][NNY][NNZ][3],
    double        UD[NNX][NNY][NNZ][3],
    double         P[NNX][NNY][NNZ],
    double        BP[NNX][NNY][NNZ][3],
    double        KX[NNX][NNY][NNZ][3]
) {
    #pragma acc kernels loop independent collapse(3) present(F, U, UD, P, BP, KX) copyin(B)
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
                    double dp1, dp2, dp3;
                    double pc0;
                    double pe1, pn1, pt1;
                    double pw1, ps1, pb1;
                    double ref, dis;

                    kx1 = KX[i    ][j    ][k    ][_X];
                    kx2 = KX[i    ][j    ][k    ][_Y];
                    kx3 = KX[i    ][j    ][k    ][_Z];
                    pc0 =  P[i    ][j    ][k    ];
                    pe1 =  P[i + 1][j    ][k    ];
                    pn1 =  P[i    ][j + 1][k    ];
                    pt1 =  P[i    ][j    ][k + 1];
                    pw1 =  P[i - 1][j    ][k    ];
                    ps1 =  P[i    ][j - 1][k    ];
                    pb1 =  P[i    ][j    ][k - 1];
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
                    b12 = f_see(B[f12], _B_P, MASK2);
                    b13 = f_see(B[f13], _B_P, MASK2);
                    b22 = f_see(B[f22], _B_P, MASK2);
                    b23 = f_see(B[f23], _B_P, MASK2);
                    b32 = f_see(B[f32], _B_P, MASK2);
                    b33 = f_see(B[f33], _B_P, MASK2);
                    if (f13) {
                        _pre_bc_eva(m13, pc0, pe1, ref, dis);
                        pe1 = 2 * bc_eva(b13, ref, dis, BP[i][j][k][_E]) - pc0;
                    }
                    if (f23) {
                        _pre_bc_eva(m23, pc0, pn1, ref, dis);
                        pn1 = 2 * bc_eva(b23, ref, dis, BP[i][j][k][_N]) - pc0;
                    }
                    if (f33) {
                        _pre_bc_eva(m33, pc0, pt1, ref, dis);
                        pt1 = 2 * bc_eva(b33, ref, dis, BP[i][j][k][_T]) - pc0;
                    }
                    if (f12) {
                        _pre_bc_eva(m12, pw1, pc0, ref, dis);
                        pw1 = 2 * bc_eva(b12, ref, dis, BP[i - 1][j][k][_E]) - pc0;
                    }
                    if (f22) {
                        _pre_bc_eva(m22, ps1, pc0, ref, dis);
                        ps1 = 2 * bc_eva(b22, ref, dis, BP[i][j - 1][k][_N]) - pc0;
                    }
                    if (f32) {
                        _pre_bc_eva(m32, pb1, pc0, ref, dis);
                        pb1 = 2 * bc_eva(b32, ref, dis, BP[i][j][k - 1][_T]) - pc0;
                    }

                    dp1 = 0.5 * (pe1 - pw1) * kx1;
                    dp2 = 0.5 * (pn1 - ps1) * kx2;
                    dp3 = 0.5 * (pt1 - pb1) * kx3;

                    U[i][j][k][_U] = UD[i][j][k][_U] - DT * dp1;
                    U[i][j][k][_V] = UD[i][j][k][_V] - DT * dp2;
                    U[i][j][k][_W] = UD[i][j][k][_W] - DT * dp3;
                }
            }
        }
    }
}

void ns_correction_f(
    unsigned int   F[NNX][NNY][NNZ],
    double         U[NNX][NNY][NNZ][3],
    double        UU[NNX][NNY][NNZ][3],
    double       UUD[NNX][NNY][NNZ][3],
    double        BU[NNX][NNY][NNZ][3][3],
    double         P[NNX][NNY][NNZ],
    double        BP[NNX][NNY][NNZ][3],
    double         X[NNX][NNY][NNZ][3],
    double        KX[NNX][NNY][NNZ][3],
    double         G[NNX][NNY][NNZ][3]
) {
    #pragma acc kernels loop independent collapse(3) present(F, U, UU, UUD, BU, P, BP, X, KX, G) copyin(B)
    for (int i = I0 - 1; i <= I1; i ++) {
        for (int j = J0 - 1; j <= J1; j ++) {
            for (int k = K0 - 1; k <= K1; k ++) {
                unsigned int ff;
                unsigned int f3, b3, m3;
                double p0, p1, dp;
                double uf;
                double x1, x2, x3;
                double k1, k2, k3;
                double g0, g1;
                double de;
                double rf, di;
                ff = F[i][j][k];

                m3 = f_see(ff, _M_E, MASK1);
                f3 = f_see(ff, _F_E, MASK8);
                b3 = f_see(B[f3], _B_U, MASK2);
                if (f3) {
                    rf = U[i + m3][j][k][_U];
                    di = 0.5 - m3;
                    uf = bc_eva(b3, rf, di, BU[i][j][k][_E][_U]);
                    x1 =         X[i + 1][j][k][_X] -  X[i][j][k][_X];
                    k2 = 0.5 * (KX[i + 1][j][k][_Y] + KX[i][j][k][_Y]);
                    k3 = 0.5 * (KX[i + 1][j][k][_Z] + KX[i][j][k][_Z]);
                    k1 = 1 / x1;
                    x2 = 1 / k2;
                    x3 = 1 / k3;
                    de = x1 * x2 * x3;
                    uf = de * k1 * uf;
                } else {
                    p0 = P[i    ][j][k];
                    p1 = P[i + 1][j][k];
                    g0 = G[i    ][j][k][_X];
                    g1 = G[i + 1][j][k][_X];
                    dp = 0.5 * (g0 + g1) * (p1 - p0);
                    uf = UUD[i][j][k][_U] - DT * dp;
                }
                UU[i][j][k][_U] = uf;

                m3 = f_see(ff, _M_N, MASK1);
                f3 = f_see(ff, _F_N, MASK8);
                b3 = f_see(B[f3], _B_V, MASK2);
                if (f3) {
                    rf = U[i][j + m3][k][_V];
                    di = 0.5 - m3;
                    uf = bc_eva(b3, rf, di, BU[i][j][k][_N][_V]);
                    k1 = 0.5 * (KX[i][j + 1][k][_X] + KX[i][j][k][_X]);
                    x2 =         X[i][j + 1][k][_Y] -  X[i][j][k][_Y];
                    k3 = 0.5 * (KX[i][j + 1][k][_Z] + KX[i][j][k][_Z]);
                    x1 = 1 / k1;
                    k2 = 1 / x2;
                    x3 = 1 / k3;
                    de = x1 * x2 * x3;
                    uf = de * k2 * uf;
                } else {
                    p0 = P[i][j    ][k];
                    p1 = P[i][j + 1][k];
                    g0 = G[i][j    ][k][_Y];
                    g1 = G[i][j + 1][k][_Y];
                    dp = 0.5 * (g0 + g1) * (p1 - p0);
                    uf = UUD[i][j][k][_V] - DT * dp;
                }
                UU[i][j][k][_V] = uf;

                m3 = f_see(ff, _M_T, MASK1);
                f3 = f_see(ff, _F_T, MASK8);
                b3 = f_see(B[f3], _B_W, MASK2);
                if (f3) {
                    rf = U[i][j][k + m3][_W];
                    di = 0.5 - m3;
                    uf = bc_eva(b3, rf, di, BU[i][j][k][_T][_W]);
                    k1 = 0.5 * (KX[i][j][k + 1][_X] + KX[i][j][k][_X]);
                    k2 = 0.5 * (KX[i][j][k + 1][_Y] + KX[i][j][k][_Y]);
                    x3 =         X[i][j][k + 1][_Z] -  X[i][j][k][_Z];
                    x1 = 1 / k1;
                    x2 = 1 / k2;
                    k3 = 1 / x3;
                    de = x1 * x2 * x3;
                    uf = de * k3 * uf;
                } else {
                    p0 = P[i][j][k    ];
                    p1 = P[i][j][k + 1];
                    g0 = G[i][j][k    ][_Z];
                    g1 = G[i][j][k + 1][_Z];
                    dp = 0.5 * (g0 + g1) * (p1 - p0);
                    uf = UUD[i][j][k][_W] - DT * dp;
                }
                UU[i][j][k][_W] = uf;
            }
        }
    }
}

void ns_correction_p(
    double  P[NNX][NNY][NNZ],
    double PP[NNX][NNY][NNZ]
) {
    #pragma acc kernels loop independent collapse(3) present(P, PP)
    for (int i = I0; i <= I1; i ++) {
        for (int j = J0; j <= J1; j ++) {
            for (int k = K0; k <= K1; k ++) {
                P[i][j][k] += PP[i][j][k];
            }
        }
    }
}
