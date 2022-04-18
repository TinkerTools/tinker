#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(kpitor, maxnpt);
extern double* TINKER_MOD(kpitor, ptcon);
extern char (*TINKER_MOD(kpitor, kpt))[8];
#ifdef __cplusplus
}
#endif
