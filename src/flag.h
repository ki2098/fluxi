#ifndef FLAG_H
#define FLAG_H

/* cell state flag : 1 bit */
#define _ACTIVE 0
/* and boundary indices flags : 8 bits */
#define _F_E    1
#define _F_N    9
#define _F_T    17
/* boundary direction flags : 1 bit */
#define _M_E    25
#define _M_N    26
#define _M_T    27
/* boundary type flags : 1 bit */
#define _D_U    0
#define _N_U    1
#define _D_V    2
#define _N_V    3
#define _D_W    4
#define _N_W    5
#define _D_P    6
#define _N_P    7
#define _BT_U   0
#define _BT_V   2
#define _BT_W   4
#define _BT_P   6
#define _BT_D   0
#define _BT_N   1

/* flag operations */
#define MASK8  255u
#define MASK2  3u
#define MASK1  1u

static unsigned int f_see(unsigned int flag, unsigned int position, unsigned int mask) {
    return (flag >> position) & mask;
}

static unsigned int f_set(unsigned int flag, unsigned int position, unsigned int value, unsigned int mask) {
    flag = flag & ~(mask << position);
    flag = flag | (value << position);
    return flag;
}

#endif
