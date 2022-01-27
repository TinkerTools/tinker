#pragma once

#include "macro.hh"
#include "sizes.hh"

namespace tinker { namespace potfit {
using namespace sizes;

extern int& nconf;
extern int& namax;
extern int& ngatm;
extern int& nfatm;
extern int (&npgrid)[maxref];
extern int*& ipgrid;
extern double& wresp;
extern double (&xdpl0)[maxref];
extern double (&ydpl0)[maxref];
extern double (&zdpl0)[maxref];
extern double (&xxqpl0)[maxref];
extern double (&xyqpl0)[maxref];
extern double (&xzqpl0)[maxref];
extern double (&yyqpl0)[maxref];
extern double (&yzqpl0)[maxref];
extern double (&zzqpl0)[maxref];
extern double*& fit0;
extern double*& fchg;
extern double*& fpol;
extern double*& fcpen;
extern double*& pgrid;
extern double*& epot;
extern int& use_dpl;
extern int& use_qpl;
extern int& fit_mpl;
extern int& fit_dpl;
extern int& fit_qpl;
extern int& fit_chgpen;
extern int*& fitchg;
extern int*& fitpol;
extern int*& fitcpen;
extern int*& gatm;
extern int*& fatm;
extern char (&resptyp)[4];
extern char (*&varpot)[6];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(potfit, nconf);
extern "C" int TINKER_MOD(potfit, namax);
extern "C" int TINKER_MOD(potfit, ngatm);
extern "C" int TINKER_MOD(potfit, nfatm);
extern "C" int TINKER_MOD(potfit, npgrid)[maxref];
extern "C" int* TINKER_MOD(potfit, ipgrid);
extern "C" double TINKER_MOD(potfit, wresp);
extern "C" double TINKER_MOD(potfit, xdpl0)[maxref];
extern "C" double TINKER_MOD(potfit, ydpl0)[maxref];
extern "C" double TINKER_MOD(potfit, zdpl0)[maxref];
extern "C" double TINKER_MOD(potfit, xxqpl0)[maxref];
extern "C" double TINKER_MOD(potfit, xyqpl0)[maxref];
extern "C" double TINKER_MOD(potfit, xzqpl0)[maxref];
extern "C" double TINKER_MOD(potfit, yyqpl0)[maxref];
extern "C" double TINKER_MOD(potfit, yzqpl0)[maxref];
extern "C" double TINKER_MOD(potfit, zzqpl0)[maxref];
extern "C" double* TINKER_MOD(potfit, fit0);
extern "C" double* TINKER_MOD(potfit, fchg);
extern "C" double* TINKER_MOD(potfit, fpol);
extern "C" double* TINKER_MOD(potfit, fcpen);
extern "C" double* TINKER_MOD(potfit, pgrid);
extern "C" double* TINKER_MOD(potfit, epot);
extern "C" int TINKER_MOD(potfit, use_dpl);
extern "C" int TINKER_MOD(potfit, use_qpl);
extern "C" int TINKER_MOD(potfit, fit_mpl);
extern "C" int TINKER_MOD(potfit, fit_dpl);
extern "C" int TINKER_MOD(potfit, fit_qpl);
extern "C" int TINKER_MOD(potfit, fit_chgpen);
extern "C" int* TINKER_MOD(potfit, fitchg);
extern "C" int* TINKER_MOD(potfit, fitpol);
extern "C" int* TINKER_MOD(potfit, fitcpen);
extern "C" int* TINKER_MOD(potfit, gatm);
extern "C" int* TINKER_MOD(potfit, fatm);
extern "C" char TINKER_MOD(potfit, resptyp)[4];
extern "C" char (*TINKER_MOD(potfit, varpot))[6];

int& nconf = TINKER_MOD(potfit, nconf);
int& namax = TINKER_MOD(potfit, namax);
int& ngatm = TINKER_MOD(potfit, ngatm);
int& nfatm = TINKER_MOD(potfit, nfatm);
int (&npgrid)[maxref] = TINKER_MOD(potfit, npgrid);
int*& ipgrid = TINKER_MOD(potfit, ipgrid);
double& wresp = TINKER_MOD(potfit, wresp);
double (&xdpl0)[maxref] = TINKER_MOD(potfit, xdpl0);
double (&ydpl0)[maxref] = TINKER_MOD(potfit, ydpl0);
double (&zdpl0)[maxref] = TINKER_MOD(potfit, zdpl0);
double (&xxqpl0)[maxref] = TINKER_MOD(potfit, xxqpl0);
double (&xyqpl0)[maxref] = TINKER_MOD(potfit, xyqpl0);
double (&xzqpl0)[maxref] = TINKER_MOD(potfit, xzqpl0);
double (&yyqpl0)[maxref] = TINKER_MOD(potfit, yyqpl0);
double (&yzqpl0)[maxref] = TINKER_MOD(potfit, yzqpl0);
double (&zzqpl0)[maxref] = TINKER_MOD(potfit, zzqpl0);
double*& fit0 = TINKER_MOD(potfit, fit0);
double*& fchg = TINKER_MOD(potfit, fchg);
double*& fpol = TINKER_MOD(potfit, fpol);
double*& fcpen = TINKER_MOD(potfit, fcpen);
double*& pgrid = TINKER_MOD(potfit, pgrid);
double*& epot = TINKER_MOD(potfit, epot);
int& use_dpl = TINKER_MOD(potfit, use_dpl);
int& use_qpl = TINKER_MOD(potfit, use_qpl);
int& fit_mpl = TINKER_MOD(potfit, fit_mpl);
int& fit_dpl = TINKER_MOD(potfit, fit_dpl);
int& fit_qpl = TINKER_MOD(potfit, fit_qpl);
int& fit_chgpen = TINKER_MOD(potfit, fit_chgpen);
int*& fitchg = TINKER_MOD(potfit, fitchg);
int*& fitpol = TINKER_MOD(potfit, fitpol);
int*& fitcpen = TINKER_MOD(potfit, fitcpen);
int*& gatm = TINKER_MOD(potfit, gatm);
int*& fatm = TINKER_MOD(potfit, fatm);
char (&resptyp)[4] = TINKER_MOD(potfit, resptyp);
char (*&varpot)[6] = TINKER_MOD(potfit, varpot);
#endif
} }
