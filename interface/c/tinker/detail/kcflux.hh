#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxncfb 2000
#define TINKER_MOD__maxncfa 2000
extern double TINKER_MOD(kcflux, cflb)[TINKER_MOD__maxncfb];
extern double TINKER_MOD(kcflux, cfla)[TINKER_MOD__maxncfa][2];
extern double TINKER_MOD(kcflux, cflab)[TINKER_MOD__maxncfa][2];
extern char TINKER_MOD(kcflux, kcfb)[TINKER_MOD__maxncfb][8];
extern char TINKER_MOD(kcflux, kcfa)[TINKER_MOD__maxncfa][12];
#ifdef __cplusplus
}
#endif
