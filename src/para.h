#ifndef PARA_H
#define PARA_H

/* fluid parameter */
static const double   RE    = 40000;
static const double   RI    = 1 / RE;

/* solver parameter */
static const double   DT    = 0.004;
static const double   EPOI  = 1E-5;
static const double   EDIV0 = 1E-2;
static const double   EDIV1 = 1E-2;
static const double   OMEGA = 1.2;
static const double   ALPHA = 1.0 / 24.0;
static const int      MAXIT = 1000;
static const int      NSTEP = 50000;

/* turbulence paratemer */
static const double   C_s   = 0.1;

#endif
