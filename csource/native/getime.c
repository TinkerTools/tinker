#include <time.h>

typedef struct
    {
    double cputim;
    } chrono_t;

extern chrono_t chrono_;

double
getime(void)
    {
    return clock()/(double)CLOCKS_PER_SEC-chrono_.cputim;
    }

int
getime_(double* elapsed)
    {
    *elapsed = getime();
    return 0;
    }
