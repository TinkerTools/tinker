#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(thrmint, tibin);
extern int TINKER_MOD(thrmint, tinbcount);
extern int TINKER_MOD(thrmint, tinbin);
extern int TINKER_MOD(thrmint, tinblock);
extern int TINKER_MOD(thrmint, tinbsave);
extern int TINKER_MOD(thrmint, tinbtot);
extern int TINKER_MOD(thrmint, tinequil);
extern int TINKER_MOD(thrmint, tinstepavg);
extern int TINKER_MOD(thrmint, tiwindow);
extern int* TINKER_MOD(thrmint, tiwinend);
extern double TINKER_MOD(thrmint, tieqratio);
extern double TINKER_MOD(thrmint, tilmda);
extern double* TINKER_MOD(thrmint, tidedllist);
extern double* TINKER_MOD(thrmint, tifraclist);
extern double* TINKER_MOD(thrmint, tilmdadedl);
extern double* TINKER_MOD(thrmint, tilmdadedlstd);
extern double* TINKER_MOD(thrmint, tilmdahist);
extern double* TINKER_MOD(thrmint, tilmdalist);
extern char TINKER_MOD(thrmint, tifile)[240];
#ifdef __cplusplus
}
#endif
