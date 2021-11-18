#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern double TINKER_MOD(warp, deform);
extern double TINKER_MOD(warp, difft);
extern double TINKER_MOD(warp, diffv);
extern double TINKER_MOD(warp, diffc);
extern double* TINKER_MOD(warp, m2);
extern int TINKER_MOD(warp, use_smooth);
extern int TINKER_MOD(warp, use_dem);
extern int TINKER_MOD(warp, use_gda);
extern int TINKER_MOD(warp, use_tophat);
extern int TINKER_MOD(warp, use_stophat);
#ifdef __cplusplus
}
#endif
