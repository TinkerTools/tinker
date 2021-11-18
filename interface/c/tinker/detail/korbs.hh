#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxnpi 500
#define TINKER_MOD__maxnpi5 200
#define TINKER_MOD__maxnpi4 200
extern double TINKER_MOD(korbs, sslope)[TINKER_MOD__maxnpi];
extern double TINKER_MOD(korbs, sslope5)[TINKER_MOD__maxnpi5];
extern double TINKER_MOD(korbs, sslope4)[TINKER_MOD__maxnpi4];
extern double TINKER_MOD(korbs, tslope)[TINKER_MOD__maxnpi];
extern double TINKER_MOD(korbs, tslope5)[TINKER_MOD__maxnpi5];
extern double TINKER_MOD(korbs, tslope4)[TINKER_MOD__maxnpi4];
extern double* TINKER_MOD(korbs, electron);
extern double* TINKER_MOD(korbs, ionize);
extern double* TINKER_MOD(korbs, repulse);
extern char TINKER_MOD(korbs, kpi)[TINKER_MOD__maxnpi][8];
extern char TINKER_MOD(korbs, kpi5)[TINKER_MOD__maxnpi5][8];
extern char TINKER_MOD(korbs, kpi4)[TINKER_MOD__maxnpi4][8];
#ifdef __cplusplus
}
#endif
