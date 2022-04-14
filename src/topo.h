#ifndef TOPO_H
#define TOPO_H

/* domain */
const double X0   = 0.0;
const double Y0   = 0.0;
const double Z0   = 0.0;
const double X1   = 40.0;
const double Y1   = 20.0;
const double Z1   = 0.05;
const int    NX   = 800;
const int    NY   = 400;
const int    NZ   = 1;
const int    NNX  = NX + 4;
const int    NNY  = NY + 4;
const int    NNZ  = NZ + 4;
const int    I0   = 2;
const int    I1   = I0 - 1 + NX;
const int    J0   = 2;
const int    J1   = J0 - 1 + NY;
const int    K0   = 2;
const int    K1   = K0 - 1 + NZ;

/* square cylinder */
const double X0sq = 9.5;
const double X1sq = 10.5;
const double Y0sq = 9.5;
const double Y1sq = 10.5;
const double Z0sq = Z0;
const double Z1sq = Z1;
const int    I0sq = 192;
const int    I1sq = 211;
const int    J0sq = 192;
const int    J1sq = 211;
const int    K0sq = K0;
const int    K1sq = K1;

#endif