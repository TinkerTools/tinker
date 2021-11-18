#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxnmp 2000
extern double TINKER_MOD(kmulti, multip)[TINKER_MOD__maxnmp][13];
extern char TINKER_MOD(kmulti, mpaxis)[TINKER_MOD__maxnmp][8];
extern char TINKER_MOD(kmulti, kmp)[TINKER_MOD__maxnmp][16];
#ifdef __cplusplus
}
#endif
