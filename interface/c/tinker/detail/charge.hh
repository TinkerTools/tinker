#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(charge, nion);
extern int* TINKER_MOD(charge, iion);
extern int* TINKER_MOD(charge, jion);
extern int* TINKER_MOD(charge, kion);
extern double* TINKER_MOD(charge, pchg);
extern double* TINKER_MOD(charge, pchg0);
#ifdef __cplusplus
}
#endif
