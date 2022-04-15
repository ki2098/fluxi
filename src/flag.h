#ifndef FLAG_H
#define FLAG_H

/* cell and boundary indices flags */
#define ACT   0
#define F_E   1
#define F_N   9
#define F_T   17

/* boundary type flags */
#define D_U   1
#define N_U   2
#define D_P   3
#define N_P   4

/* flag operations */
#define FACE  255u

unsigned int see_face(unsigned int f, int face) {
    return (f >> face) & FACE;
}

unsigned int set_face(unsigned int f, int face, unsigned int value) {
    f = f & ~(FACE << face);
    f = f | (value << face);
    return f;
}

unsigned int set_bit(unsigned int f, int position, unsigned int value) {
    f = f & ~(1u << position);
    f = f | (value << position);
    return f;
}

unsigned int see_bit(unsigned int f, int position) {
    return (f >> position) & 1u;
}

#endif