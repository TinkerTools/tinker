#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(kopbnd, maxnopb);
extern double* TINKER_MOD(kopbnd, opbn);
extern char (*TINKER_MOD(kopbnd, kopb))[16];
#ifdef __cplusplus
}
#endif
