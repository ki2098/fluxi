#ifndef PARA_H
#define PARA_H

/* fluid parameters */
static const double   RE    = 1000;
static const double   RI    = 1 / RE;

/* solver parameter */
static const double   DT    = 0.005;
static const double   EPOI  = 1E-6;
static const double   EDIV  = 1E-4;
static const double   OMEGA = 1.2;
static const double   ALPHA = 1.0 / 24.0;
static const int      MAXIT = 1000;
static const int      NSTEP = 40000;

#endif