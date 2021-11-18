#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(bndstr, nbond);
extern int* TINKER_MOD(bndstr, ibnd);
extern double* TINKER_MOD(bndstr, bk);
extern double* TINKER_MOD(bndstr, bl);
#ifdef __cplusplus
}
#endif
