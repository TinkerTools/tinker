#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(ost, fli0);
extern int TINKER_MOD(ost, iost);
extern int TINKER_MOD(ost, iosthist);
extern int TINKER_MOD(ost, nflmda);
extern int TINKER_MOD(ost, nlmda);
extern int TINKER_MOD(ost, nmetahist);
extern int TINKER_MOD(ost, nosthist);
extern int TINKER_MOD(ost, ostnavg);
extern int TINKER_MOD(ost, ostnequil);
extern int TINKER_MOD(ost, ostemexp);
extern int TINKER_MOD(ost, ostepexp);
extern int TINKER_MOD(ost, ostevexp);
extern int TINKER_MOD(ost, ostinvemn);
extern int TINKER_MOD(ost, ostinvepn);
extern int TINKER_MOD(ost, ostinvevn);
extern int TINKER_MOD(ost, sizeosthist);
extern int TINKER_MOD(ost, sizemetahist);
extern int* TINKER_MOD(ost, osthist);
extern int* TINKER_MOD(ost, ostnext);
extern int* TINKER_MOD(ost, osthead);
extern double TINKER_MOD(ost, deffdl);
extern double TINKER_MOD(ost, eosttot);
extern double TINKER_MOD(ost, hbias);
extern double TINKER_MOD(ost, maxwfhist);
extern double TINKER_MOD(ost, maxwlhist);
extern double TINKER_MOD(ost, ostddgdl);
extern double TINKER_MOD(ost, ostdedl);
extern double TINKER_MOD(ost, ostdedlavg);
extern double TINKER_MOD(ost, ostdedlstd);
extern double TINKER_MOD(ost, ostdgdl);
extern double TINKER_MOD(ost, ostdt);
extern double TINKER_MOD(ost, ostelmda0);
extern double TINKER_MOD(ost, ostelmda1);
extern double TINKER_MOD(ost, osteqratio);
extern double TINKER_MOD(ost, ostfriction);
extern double TINKER_MOD(ost, ostinvemeps);
extern double TINKER_MOD(ost, ostinvepeps);
extern double TINKER_MOD(ost, ostinveveps);
extern double TINKER_MOD(ost, ostlambda);
extern double TINKER_MOD(ost, ostlambdaavg);
extern double TINKER_MOD(ost, ostlambdastd);
extern double TINKER_MOD(ost, ostmass);
extern double TINKER_MOD(ost, ostplmda0);
extern double TINKER_MOD(ost, ostplmda1);
extern double TINKER_MOD(ost, oststdev);
extern double TINKER_MOD(ost, osttheta);
extern double TINKER_MOD(ost, ostvlmda0);
extern double TINKER_MOD(ost, ostvlmda1);
extern double TINKER_MOD(ost, ostvtheta);
extern double TINKER_MOD(ost, wfhist);
extern double TINKER_MOD(ost, wflmda);
extern double TINKER_MOD(ost, wflmda2);
extern double TINKER_MOD(ost, wlhist);
extern double TINKER_MOD(ost, wlmda);
extern double TINKER_MOD(ost, wlmda2);
extern double* TINKER_MOD(ost, fkernel);
extern double* TINKER_MOD(ost, fsumkernel);
extern double* TINKER_MOD(ost, gfkernel);
extern double* TINKER_MOD(ost, metahhist);
extern double* TINKER_MOD(ost, metalhist);
extern double* TINKER_MOD(ost, metawhist);
extern double* TINKER_MOD(ost, ostflist);
extern double* TINKER_MOD(ost, ostfhist);
extern double* TINKER_MOD(ost, osthhist);
extern double* TINKER_MOD(ost, ostllist);
extern double* TINKER_MOD(ost, ostlhist);
extern double* TINKER_MOD(ost, ostwfhist);
extern double* TINKER_MOD(ost, ostwlhist);
extern double* TINKER_MOD(ost, gkernel);
extern double* TINKER_MOD(ost, glfkernel);
extern double* TINKER_MOD(ost, glkernel);
extern double* TINKER_MOD(ost, pfkernel);
extern int TINKER_MOD(ost, fastkernel);
extern int TINKER_MOD(ost, metarestart);
extern int TINKER_MOD(ost, ostinterpol);
extern int TINKER_MOD(ost, ostrestart);
extern int TINKER_MOD(ost, use_meta);
extern int TINKER_MOD(ost, use_metadyn);
extern int TINKER_MOD(ost, use_ost);
extern int TINKER_MOD(ost, use_ostdyn);
extern int TINKER_MOD(ost, use_pol4f);
extern int TINKER_MOD(ost, use_pol4i);
extern char TINKER_MOD(ost, ostemap)[3];
extern char TINKER_MOD(ost, ostpmap)[3];
extern char TINKER_MOD(ost, ostvmap)[3];
#ifdef __cplusplus
}
#endif
