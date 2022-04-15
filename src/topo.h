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

void topo_init(
    int F[NNX][NNY][NNZ]
);

#endif