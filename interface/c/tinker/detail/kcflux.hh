#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(kcflux, maxncfb);
extern int TINKER_MOD(kcflux, maxncfa);
extern double* TINKER_MOD(kcflux, cflb);
extern double* TINKER_MOD(kcflux, cfla);
extern double* TINKER_MOD(kcflux, cflab);
extern char (*TINKER_MOD(kcflux, kcfb))[8];
extern char (*TINKER_MOD(kcflux, kcfa))[12];
#ifdef __cplusplus
}
#endif
