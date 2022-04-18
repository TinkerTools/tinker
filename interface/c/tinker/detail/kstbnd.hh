#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(kstbnd, maxnsb);
extern double* TINKER_MOD(kstbnd, stbn);
extern char (*TINKER_MOD(kstbnd, ksb))[12];
#ifdef __cplusplus
}
#endif
