#ifndef FLAG_H
#define FLAG_H

#define Active 0
#define B_e    1
#define B_n    2
#define B_t    3

#define D_e    1
#define N_e    2
#define D_n    3
#define N_n    4
#define D_t    5
#define N_t    6

unsigned int set_flag(unsigned int f, int position) {
    return f | (1u << position);
}

unsigned int see_flag(unsigned int f, int position) {
    return (f >> position) & 1u;
}

unsigned int clear_flag(unsigned int f, int position) {
    return f & ~(1u << position);
}

#endif