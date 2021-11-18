#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxnopd 500
extern double TINKER_MOD(kopdst, opds)[TINKER_MOD__maxnopd];
extern char TINKER_MOD(kopdst, kopd)[TINKER_MOD__maxnopd][16];
#ifdef __cplusplus
}
#endif
