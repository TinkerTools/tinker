#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(angbnd, nangle);
extern int* TINKER_MOD(angbnd, iang);
extern double* TINKER_MOD(angbnd, ak);
extern double* TINKER_MOD(angbnd, anat);
extern double* TINKER_MOD(angbnd, afld);
#ifdef __cplusplus
}
#endif
