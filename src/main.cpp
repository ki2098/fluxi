#include <stdio.h>
#include "flag.h"

int main(void) {
    unsigned int f = 0;
    f = set_flag(f, Active);
    f = set_flag(f, N_e);
    printf("%u\n", see_flag(f, Active));
    printf("%u\n", see_flag(f, N_e));
    printf("%u\n", see_flag(f, N_n));
    f = clear_flag(f, N_e);
    printf("%u\n", see_flag(f, Active));
    printf("%u\n", see_flag(f, N_e));
    printf("%u\n", see_flag(f, N_n));

    return 0;
}