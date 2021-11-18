#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxion 10
extern int TINKER_MOD(pbstuf, ionn);
extern int TINKER_MOD(pbstuf, dime)[3];
extern int TINKER_MOD(pbstuf, ionq)[TINKER_MOD__maxion];
extern double TINKER_MOD(pbstuf, pbe);
extern double TINKER_MOD(pbstuf, pdie);
extern double TINKER_MOD(pbstuf, sdie);
extern double TINKER_MOD(pbstuf, srad);
extern double TINKER_MOD(pbstuf, swin);
extern double TINKER_MOD(pbstuf, sdens);
extern double TINKER_MOD(pbstuf, smin);
extern double TINKER_MOD(pbstuf, grid)[3];
extern double TINKER_MOD(pbstuf, gcent)[3];
extern double TINKER_MOD(pbstuf, cgrid)[3];
extern double TINKER_MOD(pbstuf, cgcent)[3];
extern double TINKER_MOD(pbstuf, fgrid)[3];
extern double TINKER_MOD(pbstuf, fgcent)[3];
extern double TINKER_MOD(pbstuf, ionr)[TINKER_MOD__maxion];
extern double TINKER_MOD(pbstuf, ionc)[TINKER_MOD__maxion];
extern double* TINKER_MOD(pbstuf, apbe);
extern double* TINKER_MOD(pbstuf, pbep);
extern double* TINKER_MOD(pbstuf, pbfp);
extern double* TINKER_MOD(pbstuf, pbtp);
extern double* TINKER_MOD(pbstuf, pbeuind);
extern double* TINKER_MOD(pbstuf, pbeuinp);
extern char TINKER_MOD(pbstuf, pbtyp)[20];
extern char TINKER_MOD(pbstuf, pbsoln)[20];
extern char TINKER_MOD(pbstuf, bcfl)[20];
extern char TINKER_MOD(pbstuf, chgm)[20];
extern char TINKER_MOD(pbstuf, srfm)[20];
#ifdef __cplusplus
}
#endif
