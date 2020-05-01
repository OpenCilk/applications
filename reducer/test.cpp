#include "opencilk_reducer.hpp"
#include <cilk/cilk.h>
#include <stdio.h>

struct Long {
    long value = 0L;

    Long() {}
    Long(long init) : value(init) {}

    void reduce(Long *right) { value += right->value; }
};

int main(int argc, char *argv[]) {
    long count = 100 * 1000;
    if (argc > 1)
        count = atol(argv[1]);
    opencilk_reducer<Long> sum(100);
    cilk_for(long i = 0; i < count; ++i)++ sum.view().value;
    cilk_sync;
    printf("sum %ld =? %ld + 100\n", sum.view().value, count);
}
