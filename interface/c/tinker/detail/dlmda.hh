#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(dlmda, elmdaexp);
extern int TINKER_MOD(dlmda, elmdainvn);
extern int TINKER_MOD(dlmda, emdtexp);
extern int TINKER_MOD(dlmda, epdtexp);
extern int TINKER_MOD(dlmda, evdtexp);
extern int TINKER_MOD(dlmda, plmdaexp);
extern int TINKER_MOD(dlmda, plmdainvn);
extern int TINKER_MOD(dlmda, vlmdaexp);
extern int TINKER_MOD(dlmda, vlmdainvn);
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
extern double TINKER_MOD(dlmda, elmdainveps);
extern double TINKER_MOD(dlmda, plmdainveps);
extern double TINKER_MOD(dlmda, qntelmda0);
extern double TINKER_MOD(dlmda, qntelmda1);
extern double TINKER_MOD(dlmda, qntplmda0);
extern double TINKER_MOD(dlmda, qntplmda1);
extern double TINKER_MOD(dlmda, qntvlmda0);
extern double TINKER_MOD(dlmda, qntvlmda1);
extern double TINKER_MOD(dlmda, vlmdainveps);
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
extern int TINKER_MOD(dlmda, use_meta);
extern int TINKER_MOD(dlmda, use_metadyn);
extern int TINKER_MOD(dlmda, use_ost);
extern int TINKER_MOD(dlmda, use_ostdyn);
extern int TINKER_MOD(dlmda, use_plmda);
extern int TINKER_MOD(dlmda, use_pol4f);
extern int TINKER_MOD(dlmda, use_pol4i);
extern int* TINKER_MOD(dlmda, douindorig);
extern char TINKER_MOD(dlmda, elmdamap)[3];
extern char TINKER_MOD(dlmda, plmdamap)[3];
extern char TINKER_MOD(dlmda, vlmdamap)[3];
#ifdef __cplusplus
}
#endif
