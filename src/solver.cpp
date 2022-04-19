#include "solver.h"
#include "bc.h"
#include <math.h>

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

void solver_sor(
    unsigned int  F[NNX][NNY][NNZ],
    double        P[NNX][NNY][NNZ],
    double       BP[NNX][NNY][NNZ][3],
    double       PS[NNX][NNY][NNZ],
    double        C[NNX][NNY][NNZ][6],
    double       &res
) {
    double err = 0;
    int    cnt = 0;

    for (int i = I0; i <= I1; i ++) {
        for (int j = J0; j <= J1; j ++) {
            for (int k = K0; k <= K1; k ++) {
                if (f_see(F[i][j][k], ACTIVE, MASK1) && (i + j + k) % 2 == 0) {
                    unsigned int f12, f13;
                    unsigned int f22, f23;
                    unsigned int f32, f33;
                    unsigned int m12, m13;
                    unsigned int m22, m23;
                    unsigned int m32, m33;
                    unsigned int b12, b13;
                    unsigned int b22, b23;
                    unsigned int b32, b33;
                    double ac0;
                    double ae1, an1, at1;
                    double aw1, as1, ab1;
                    double pc0;
                    double pe1, pn1, pt1;
                    double pw1, ps1, pb1;
                    double c01, c02, c03;
                    double c07, c08, c09;
                    double ref, dis;
                    double rhs, rc0;

                    rhs = PS[i    ][j    ][k    ] / DT;
                    c01 =  C[i    ][j    ][k    ][0];
                    c02 =  C[i    ][j    ][k    ][1];
                    c03 =  C[i    ][j    ][k    ][2];
                    c07 =  C[i    ][j    ][k    ][3];
                    c08 =  C[i    ][j    ][k    ][4];
                    c09 =  C[i    ][j    ][k    ][5];
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
                    m12 = f_see(f12, M_E, MASK1);
                    m13 = f_see(f13, M_E, MASK1);
                    m22 = f_see(f22, M_N, MASK1);
                    m23 = f_see(f23, M_N, MASK1);
                    m32 = f_see(f32, M_T, MASK1);
                    m33 = f_see(f33, M_T, MASK1);
                    f12 = f_see(f12, F_E, MASK8);
                    f13 = f_see(f13, F_E, MASK8);
                    f22 = f_see(f22, F_N, MASK8);
                    f23 = f_see(f23, F_N, MASK8);
                    f32 = f_see(f32, F_T, MASK8);
                    f33 = f_see(f33, F_T, MASK8);
                    b12 = f_see(B[f12], BT_P, MASK2);
                    b13 = f_see(B[f13], BT_P, MASK2);
                    b22 = f_see(B[f22], BT_P, MASK2);
                    b23 = f_see(B[f23], BT_P, MASK2);
                    b32 = f_see(B[f32], BT_P, MASK2);
                    b33 = f_see(B[f33], BT_P, MASK2);
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

                    ac0 = - 2 * (c01 + c02 + c03);
                    ae1 = c01 - 0.5 * c07;
                    an1 = c02 - 0.5 * c08;
                    at1 = c03 - 0.5 * c09;
                    aw1 = c01 + 0.5 * c07;
                    as1 = c02 + 0.5 * c08;
                    ab1 = c03 + 0.5 * c09;
                    rc0 = (rhs - ae1 * pe1 - an1 * pn1 - at1 * pt1 - aw1 * pw1 - as1 * ps1 - ab1 * pb1) / ac0 - pc0;
                    err = err + rc0 * rc0;

                    P[i][j][k] = pc0 + OMEGA * rc0;
                    cnt += 1;
                }
            }
        }
    }

    for (int i = I0; i <= I1; i ++) {
        for (int j = J0; j <= J1; j ++) {
            for (int k = K0; k <= K1; k ++) {
                if (f_see(F[i][j][k], ACTIVE, MASK1) && (i + j + k) % 2 == 1) {
                    unsigned int f12, f13;
                    unsigned int f22, f23;
                    unsigned int f32, f33;
                    unsigned int m12, m13;
                    unsigned int m22, m23;
                    unsigned int m32, m33;
                    unsigned int b12, b13;
                    unsigned int b22, b23;
                    unsigned int b32, b33;
                    double ac0;
                    double ae1, an1, at1;
                    double aw1, as1, ab1;
                    double pc0;
                    double pe1, pn1, pt1;
                    double pw1, ps1, pb1;
                    double c01, c02, c03;
                    double c07, c08, c09;
                    double ref, dis;
                    double rhs, rc0;

                    rhs = PS[i    ][j    ][k    ] / DT;
                    c01 =  C[i    ][j    ][k    ][0];
                    c02 =  C[i    ][j    ][k    ][1];
                    c03 =  C[i    ][j    ][k    ][2];
                    c07 =  C[i    ][j    ][k    ][3];
                    c08 =  C[i    ][j    ][k    ][4];
                    c09 =  C[i    ][j    ][k    ][5];
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
                    m12 = f_see(f12, M_E, MASK1);
                    m13 = f_see(f13, M_E, MASK1);
                    m22 = f_see(f22, M_N, MASK1);
                    m23 = f_see(f23, M_N, MASK1);
                    m32 = f_see(f32, M_T, MASK1);
                    m33 = f_see(f33, M_T, MASK1);
                    f12 = f_see(f12, F_E, MASK8);
                    f13 = f_see(f13, F_E, MASK8);
                    f22 = f_see(f22, F_N, MASK8);
                    f23 = f_see(f23, F_N, MASK8);
                    f32 = f_see(f32, F_T, MASK8);
                    f33 = f_see(f33, F_T, MASK8);
                    b12 = f_see(B[f12], BT_P, MASK2);
                    b13 = f_see(B[f13], BT_P, MASK2);
                    b22 = f_see(B[f22], BT_P, MASK2);
                    b23 = f_see(B[f23], BT_P, MASK2);
                    b32 = f_see(B[f32], BT_P, MASK2);
                    b33 = f_see(B[f33], BT_P, MASK2);
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

                    ac0 = - 2 * (c01 + c02 + c03);
                    ae1 = c01 - 0.5 * c07;
                    an1 = c02 - 0.5 * c08;
                    at1 = c03 - 0.5 * c09;
                    aw1 = c01 + 0.5 * c07;
                    as1 = c02 + 0.5 * c08;
                    ab1 = c03 + 0.5 * c09;
                    rc0 = (rhs - ae1 * pe1 - an1 * pn1 - at1 * pt1 - aw1 * pw1 - as1 * ps1 - ab1 * pb1) / ac0 - pc0;
                    err = err + rc0 * rc0;

                    P[i][j][k] = pc0 + OMEGA * rc0;
                    cnt += 1;
                }
            }
        }
    }

    res = sqrt(err / cnt);
}