#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern double TINKER_MOD(paths, pnorm);
extern double TINKER_MOD(paths, acoeff)[7][7];
extern double* TINKER_MOD(paths, pc0);
extern double* TINKER_MOD(paths, pc1);
extern double* TINKER_MOD(paths, pvect);
extern double* TINKER_MOD(paths, pstep);
extern double* TINKER_MOD(paths, pzet);
extern double* TINKER_MOD(paths, gc);
#ifdef __cplusplus
}
#endif
