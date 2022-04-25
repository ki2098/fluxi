#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "flag.h"

/* number of boundaries */
#define NB      4 
#define NOTHING 0
#define WALL    1
#define INFLOW  2
#define OUTFLOW 3
#define SLIP    4

/* predefined boundary values */
#define UINFLOW 1.0
#define VINFLOW 0.0
#define WINFLOW 0.0

/* boundary definitions */
static const unsigned int B[NB + 1] = {
/*  u                          v                          w                          p */
    0                        | 0                        | 0                        | 0                       , // no boundary
    f_set(0, _D_U, 1, MASK1) | f_set(0, _D_V, 1, MASK1) | f_set(0, _D_W, 1, MASK1) | f_set(0, _N_P, 1, MASK1), // wall
    f_set(0, _D_U, 1, MASK1) | f_set(0, _D_V, 1, MASK1) | f_set(0, _D_W, 1, MASK1) | f_set(0, _N_P, 1, MASK1), // inflow
    f_set(0, _D_U, 1, MASK1) | f_set(0, _D_V, 1, MASK1) | f_set(0, _D_W, 1, MASK1) | f_set(0, _N_P, 1, MASK1), // outflow
    f_set(0, _N_U, 1, MASK1) | f_set(0, _D_V, 1, MASK1) | f_set(0, _N_W, 1, MASK1) | f_set(0, _N_P, 1, MASK1)  // slip (normal in y)
};


#endif
