#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(thrmint, tinbcount);
extern int TINKER_MOD(thrmint, tinbsave);
extern int TINKER_MOD(thrmint, tinbtot);
extern int TINKER_MOD(thrmint, tinstepavg);
extern double* TINKER_MOD(thrmint, tidedllist);
extern double* TINKER_MOD(thrmint, tilmdadedl);
extern double* TINKER_MOD(thrmint, tilmdadedlstd);
extern double* TINKER_MOD(thrmint, tilmdahist);
#ifdef __cplusplus
}
#endif
