#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(expol, nexpol);
extern double* TINKER_MOD(expol, kpep);
extern double* TINKER_MOD(expol, prepep);
extern double* TINKER_MOD(expol, dmppep);
extern double* TINKER_MOD(expol, polscale);
extern double* TINKER_MOD(expol, invpolscale);
extern int* TINKER_MOD(expol, lpep);
#ifdef __cplusplus
}
#endif
