#ifndef BOUNDARY2_H
#define BOUNDARY2_H

#include "flag.h"

/* number of boundaries */
#define NB      5 
#define NOTHING 0
#define WALL    1
#define INFLOW  2
#define OUTFLOW 3
#define SLIP    4

const int UINFLOW = 1.0;
const int VINFLOW = 0.0;
const int WINFLOW = 0.0;

/* boundary definitions */
const unsigned int B[NB] = {
/*  u                    p               */
    0                  | 0                 , // no boundary
    set_bit(0, D_U, 1) | set_bit(0, N_P, 1), // wall
    set_bit(0, D_U, 1) | set_bit(0, N_P, 1), // inflow
    set_bit(0, D_U, 1) | set_bit(0, N_P, 1), // outflow
    set_bit(0, D_U, 1) | set_bit(0, N_P, 1), // slip
};


#endif