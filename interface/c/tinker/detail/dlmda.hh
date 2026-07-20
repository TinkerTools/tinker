#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(dlmda, emdtexp);
extern int TINKER_MOD(dlmda, epdtexp);
extern int TINKER_MOD(dlmda, evdtexp);
extern double TINKER_MOD(dlmda, d2edl2);
extern double TINKER_MOD(dlmda, d2eldlmda2);
extern double TINKER_MOD(dlmda, d2emdl2);
extern double TINKER_MOD(dlmda, d2epdl2);
extern double TINKER_MOD(dlmda, d2evdl2);
extern double TINKER_MOD(dlmda, d2pldlmda2);
extern double TINKER_MOD(dlmda, d2vldlmda2);
extern double TINKER_MOD(dlmda, dedl);
extern double TINKER_MOD(dlmda, deldlmda);
extern double TINKER_MOD(dlmda, demdl);
extern double TINKER_MOD(dlmda, depdl);
extern double TINKER_MOD(dlmda, devdl);
extern double TINKER_MOD(dlmda, dpldlmda);
extern double TINKER_MOD(dlmda, dvldlmda);
extern double TINKER_MOD(dlmda, plambda);
extern double TINKER_MOD(dlmda, demvirdl)[3][3];
extern double TINKER_MOD(dlmda, depvirdl)[3][3];
extern double TINKER_MOD(dlmda, devvirdl)[3][3];
extern double TINKER_MOD(dlmda, dvirdl)[3][3];
extern double* TINKER_MOD(dlmda, abflxorig);
extern double* TINKER_MOD(dlmda, aflxorig);
extern double* TINKER_MOD(dlmda, bdplorig);
extern double* TINKER_MOD(dlmda, bflxorig);
extern double* TINKER_MOD(dlmda, dfmdl);
extern double* TINKER_MOD(dlmda, dfpdl);
extern double* TINKER_MOD(dlmda, dfsumdl);
extern double* TINKER_MOD(dlmda, dfvdl);
extern double* TINKER_MOD(dlmda, lcmp);
extern double* TINKER_MOD(dlmda, lcphi);
extern double* TINKER_MOD(dlmda, lfmp);
extern double* TINKER_MOD(dlmda, lfphi);
extern double* TINKER_MOD(dlmda, lqgrid);
extern double* TINKER_MOD(dlmda, pchg0orig);
extern double* TINKER_MOD(dlmda, pchgorig);
extern double* TINKER_MOD(dlmda, pcoreorig);
extern double* TINKER_MOD(dlmda, polarityorig);
extern double* TINKER_MOD(dlmda, poleorig);
extern double* TINKER_MOD(dlmda, pval0orig);
extern double* TINKER_MOD(dlmda, pvalorig);
extern int TINKER_MOD(dlmda, use_dlmda);
extern int TINKER_MOD(dlmda, use_emdt);
extern int TINKER_MOD(dlmda, use_epdt);
extern int TINKER_MOD(dlmda, use_evdt);
extern int TINKER_MOD(dlmda, use_plmda);
extern int* TINKER_MOD(dlmda, douindorig);
#ifdef __cplusplus
}
#endif
