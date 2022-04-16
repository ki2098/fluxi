#include <stdio.h>
#include "tp.h"
#include "var.h"

int main(void) {
    tp_f_init(F);
    tp_x_init(X, KX, J, G, C);

    printf("%u\n", B[WALL]);


    return 0;
}