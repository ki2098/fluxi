#ifndef FLAG_H
#define FLAG_H

/* cell state flag : 1 bit */
#define ACTIVE 0
/* and boundary indices flags : 8 bits */
#define F_E    1
#define F_N    9
#define F_T    17
/* boundary direction flags : 1 bit */
#define M_E    25
#define M_N    26
#define M_T    27
/* boundary type flags : 1 bit */
#define D_U    1
#define N_U    2
#define D_V    3
#define N_V    4
#define D_W    5
#define N_W    6
#define D_P    7
#define N_P    8
#define BT_U   1
#define BT_V   3
#define BT_W   5
#define BT_P   7
#define BT_D   0
#define BT_N   1

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