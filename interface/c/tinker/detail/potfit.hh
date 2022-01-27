#pragma once

#include "macro.hh"
#include "sizes.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(potfit, nconf);
extern int TINKER_MOD(potfit, namax);
extern int TINKER_MOD(potfit, ngatm);
extern int TINKER_MOD(potfit, nfatm);
extern int TINKER_MOD(potfit, npgrid)[TINKER_MOD__maxref];
extern int* TINKER_MOD(potfit, ipgrid);
extern double TINKER_MOD(potfit, wresp);
extern double TINKER_MOD(potfit, xdpl0)[TINKER_MOD__maxref];
extern double TINKER_MOD(potfit, ydpl0)[TINKER_MOD__maxref];
extern double TINKER_MOD(potfit, zdpl0)[TINKER_MOD__maxref];
extern double TINKER_MOD(potfit, xxqpl0)[TINKER_MOD__maxref];
extern double TINKER_MOD(potfit, xyqpl0)[TINKER_MOD__maxref];
extern double TINKER_MOD(potfit, xzqpl0)[TINKER_MOD__maxref];
extern double TINKER_MOD(potfit, yyqpl0)[TINKER_MOD__maxref];
extern double TINKER_MOD(potfit, yzqpl0)[TINKER_MOD__maxref];
extern double TINKER_MOD(potfit, zzqpl0)[TINKER_MOD__maxref];
extern double* TINKER_MOD(potfit, fit0);
extern double* TINKER_MOD(potfit, fchg);
extern double* TINKER_MOD(potfit, fpol);
extern double* TINKER_MOD(potfit, fcpen);
extern double* TINKER_MOD(potfit, pgrid);
extern double* TINKER_MOD(potfit, epot);
extern int TINKER_MOD(potfit, use_dpl);
extern int TINKER_MOD(potfit, use_qpl);
extern int TINKER_MOD(potfit, fit_mpl);
extern int TINKER_MOD(potfit, fit_dpl);
extern int TINKER_MOD(potfit, fit_qpl);
extern int TINKER_MOD(potfit, fit_chgpen);
extern int* TINKER_MOD(potfit, fitchg);
extern int* TINKER_MOD(potfit, fitpol);
extern int* TINKER_MOD(potfit, fitcpen);
extern int* TINKER_MOD(potfit, gatm);
extern int* TINKER_MOD(potfit, fatm);
extern char TINKER_MOD(potfit, resptyp)[4];
extern char (*TINKER_MOD(potfit, varpot))[6];
#ifdef __cplusplus
}
#endif
