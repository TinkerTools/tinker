#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int* TINKER_MOD(polpcg, mindex);
extern double TINKER_MOD(polpcg, pcgpeek);
extern double* TINKER_MOD(polpcg, minv);
extern int TINKER_MOD(polpcg, pcgprec);
extern int TINKER_MOD(polpcg, pcgguess);
#ifdef __cplusplus
}
#endif
