#pragma once

#include "macro.hh"
#include "sizes.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(zcoord, iz)[TINKER_MOD__maxatm][4];
extern double TINKER_MOD(zcoord, zbond)[TINKER_MOD__maxatm];
extern double TINKER_MOD(zcoord, zang)[TINKER_MOD__maxatm];
extern double TINKER_MOD(zcoord, ztors)[TINKER_MOD__maxatm];
#ifdef __cplusplus
}
#endif
