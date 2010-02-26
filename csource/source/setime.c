#include <time.h>

typedef struct
    {
    double cputim;
    } chrono_t;

extern chrono_t chrono_;
chrono_t chrono_;

void
setime(void)
    {
    chrono_.cputim=clock()/(double)CLOCKS_PER_SEC;
    }

int
setime_(void)
    {
    setime();
    return 0;
    }
