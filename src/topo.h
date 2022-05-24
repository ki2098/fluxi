#ifndef TOPO_H
#define TOPO_H

#define _E 0
#define _N 1
#define _T 2
#define _X 0
#define _Y 1
#define _Z 2
#define _U 0
#define _V 1
#define _W 2

/* domain */
static const double X0   = 0.0;
static const double Y0   = 0.0;
static const double Z0   = 0.0;
static const double X1   = 14.0;
static const double Y1   = 6.4;
static const double Z1   = 3.0;
static const int    NX   = 200;
static const int    NY   = 150;
static const int    NZ   = 120;
static const int    NNX  = NX + 4;
static const int    NNY  = NY + 4;
static const int    NNZ  = NZ + 4;
static const int    I0   = 2;
static const int    I1   = I0 - 1 + NX;
static const int    J0   = 2;
static const int    J1   = J0 - 1 + NY;
static const int    K0   = 2;
static const int    K1   = K0 - 1 + NZ;

/* square cylinder */
static const double X0sq = 3.0;
static const double X1sq = 4.0;
static const double Y0sq = 2.7;
static const double Y1sq = 3.7;
static const double Z0sq = Z0;
static const double Z1sq = 1.0;
static const int    I0sq = 52;
static const int    I1sq = 101;
static const int    J0sq = 52;
static const int    J1sq = 101;
static const int    K0sq = K0;
static const int    K1sq = 71;

#endif
