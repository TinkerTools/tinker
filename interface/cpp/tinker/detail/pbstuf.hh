#pragma once

#include "macro.hh"

namespace tinker { namespace pbstuf {
const int maxion = 10;
extern int& ionn;
extern int (&dime)[3];
extern int (&ionq)[maxion];
extern double& pbe;
extern double& pdie;
extern double& sdie;
extern double& srad;
extern double& swin;
extern double& sdens;
extern double& smin;
extern double (&grid)[3];
extern double (&gcent)[3];
extern double (&cgrid)[3];
extern double (&cgcent)[3];
extern double (&fgrid)[3];
extern double (&fgcent)[3];
extern double (&ionr)[maxion];
extern double (&ionc)[maxion];
extern double*& apbe;
extern double*& pbep;
extern double*& pbfp;
extern double*& pbtp;
extern double*& pbeuind;
extern double*& pbeuinp;
extern char (&pbtyp)[20];
extern char (&pbsoln)[20];
extern char (&bcfl)[20];
extern char (&chgm)[20];
extern char (&srfm)[20];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(pbstuf, ionn);
extern "C" int TINKER_MOD(pbstuf, dime)[3];
extern "C" int TINKER_MOD(pbstuf, ionq)[maxion];
extern "C" double TINKER_MOD(pbstuf, pbe);
extern "C" double TINKER_MOD(pbstuf, pdie);
extern "C" double TINKER_MOD(pbstuf, sdie);
extern "C" double TINKER_MOD(pbstuf, srad);
extern "C" double TINKER_MOD(pbstuf, swin);
extern "C" double TINKER_MOD(pbstuf, sdens);
extern "C" double TINKER_MOD(pbstuf, smin);
extern "C" double TINKER_MOD(pbstuf, grid)[3];
extern "C" double TINKER_MOD(pbstuf, gcent)[3];
extern "C" double TINKER_MOD(pbstuf, cgrid)[3];
extern "C" double TINKER_MOD(pbstuf, cgcent)[3];
extern "C" double TINKER_MOD(pbstuf, fgrid)[3];
extern "C" double TINKER_MOD(pbstuf, fgcent)[3];
extern "C" double TINKER_MOD(pbstuf, ionr)[maxion];
extern "C" double TINKER_MOD(pbstuf, ionc)[maxion];
extern "C" double* TINKER_MOD(pbstuf, apbe);
extern "C" double* TINKER_MOD(pbstuf, pbep);
extern "C" double* TINKER_MOD(pbstuf, pbfp);
extern "C" double* TINKER_MOD(pbstuf, pbtp);
extern "C" double* TINKER_MOD(pbstuf, pbeuind);
extern "C" double* TINKER_MOD(pbstuf, pbeuinp);
extern "C" char TINKER_MOD(pbstuf, pbtyp)[20];
extern "C" char TINKER_MOD(pbstuf, pbsoln)[20];
extern "C" char TINKER_MOD(pbstuf, bcfl)[20];
extern "C" char TINKER_MOD(pbstuf, chgm)[20];
extern "C" char TINKER_MOD(pbstuf, srfm)[20];

int& ionn = TINKER_MOD(pbstuf, ionn);
int (&dime)[3] = TINKER_MOD(pbstuf, dime);
int (&ionq)[maxion] = TINKER_MOD(pbstuf, ionq);
double& pbe = TINKER_MOD(pbstuf, pbe);
double& pdie = TINKER_MOD(pbstuf, pdie);
double& sdie = TINKER_MOD(pbstuf, sdie);
double& srad = TINKER_MOD(pbstuf, srad);
double& swin = TINKER_MOD(pbstuf, swin);
double& sdens = TINKER_MOD(pbstuf, sdens);
double& smin = TINKER_MOD(pbstuf, smin);
double (&grid)[3] = TINKER_MOD(pbstuf, grid);
double (&gcent)[3] = TINKER_MOD(pbstuf, gcent);
double (&cgrid)[3] = TINKER_MOD(pbstuf, cgrid);
double (&cgcent)[3] = TINKER_MOD(pbstuf, cgcent);
double (&fgrid)[3] = TINKER_MOD(pbstuf, fgrid);
double (&fgcent)[3] = TINKER_MOD(pbstuf, fgcent);
double (&ionr)[maxion] = TINKER_MOD(pbstuf, ionr);
double (&ionc)[maxion] = TINKER_MOD(pbstuf, ionc);
double*& apbe = TINKER_MOD(pbstuf, apbe);
double*& pbep = TINKER_MOD(pbstuf, pbep);
double*& pbfp = TINKER_MOD(pbstuf, pbfp);
double*& pbtp = TINKER_MOD(pbstuf, pbtp);
double*& pbeuind = TINKER_MOD(pbstuf, pbeuind);
double*& pbeuinp = TINKER_MOD(pbstuf, pbeuinp);
char (&pbtyp)[20] = TINKER_MOD(pbstuf, pbtyp);
char (&pbsoln)[20] = TINKER_MOD(pbstuf, pbsoln);
char (&bcfl)[20] = TINKER_MOD(pbstuf, bcfl);
char (&chgm)[20] = TINKER_MOD(pbstuf, chgm);
char (&srfm)[20] = TINKER_MOD(pbstuf, srfm);
#endif
} }
