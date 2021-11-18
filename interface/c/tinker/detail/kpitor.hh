#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxnpt 500
extern double TINKER_MOD(kpitor, ptcon)[TINKER_MOD__maxnpt];
extern char TINKER_MOD(kpitor, kpt)[TINKER_MOD__maxnpt][8];
#ifdef __cplusplus
}
#endif
