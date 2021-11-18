#pragma once

#include "macro.hh"
#include "sizes.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(refer, nref)[TINKER_MOD__maxref];
extern int TINKER_MOD(refer, refltitle)[TINKER_MOD__maxref];
extern int TINKER_MOD(refer, refleng)[TINKER_MOD__maxref];
extern int* TINKER_MOD(refer, reftyp);
extern int* TINKER_MOD(refer, n12ref);
extern int* TINKER_MOD(refer, i12ref);
extern double TINKER_MOD(refer, xboxref)[TINKER_MOD__maxref];
extern double TINKER_MOD(refer, yboxref)[TINKER_MOD__maxref];
extern double TINKER_MOD(refer, zboxref)[TINKER_MOD__maxref];
extern double TINKER_MOD(refer, alpharef)[TINKER_MOD__maxref];
extern double TINKER_MOD(refer, betaref)[TINKER_MOD__maxref];
extern double TINKER_MOD(refer, gammaref)[TINKER_MOD__maxref];
extern double* TINKER_MOD(refer, xref);
extern double* TINKER_MOD(refer, yref);
extern double* TINKER_MOD(refer, zref);
extern char (*TINKER_MOD(refer, refnam))[3];
extern char TINKER_MOD(refer, reffile)[TINKER_MOD__maxref][240];
extern char TINKER_MOD(refer, reftitle)[TINKER_MOD__maxref][240];
#ifdef __cplusplus
}
#endif
