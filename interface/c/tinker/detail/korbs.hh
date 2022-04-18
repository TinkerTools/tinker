#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(korbs, maxnpi);
extern int TINKER_MOD(korbs, maxnpi5);
extern int TINKER_MOD(korbs, maxnpi4);
extern double* TINKER_MOD(korbs, electron);
extern double* TINKER_MOD(korbs, ionize);
extern double* TINKER_MOD(korbs, repulse);
extern double* TINKER_MOD(korbs, sslope);
extern double* TINKER_MOD(korbs, sslope5);
extern double* TINKER_MOD(korbs, sslope4);
extern double* TINKER_MOD(korbs, tslope);
extern double* TINKER_MOD(korbs, tslope5);
extern double* TINKER_MOD(korbs, tslope4);
extern char (*TINKER_MOD(korbs, kpi))[8];
extern char (*TINKER_MOD(korbs, kpi5))[8];
extern char (*TINKER_MOD(korbs, kpi4))[8];
#ifdef __cplusplus
}
#endif
