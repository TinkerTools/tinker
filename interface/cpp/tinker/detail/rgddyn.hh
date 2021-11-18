#pragma once

#include "macro.hh"

namespace tinker { namespace rgddyn {
extern double*& xcmo;
extern double*& ycmo;
extern double*& zcmo;
extern double*& vcm;
extern double*& wcm;
extern double*& lm;
extern double*& vc;
extern double*& wc;
extern int*& linear;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(rgddyn, xcmo);
extern "C" double* TINKER_MOD(rgddyn, ycmo);
extern "C" double* TINKER_MOD(rgddyn, zcmo);
extern "C" double* TINKER_MOD(rgddyn, vcm);
extern "C" double* TINKER_MOD(rgddyn, wcm);
extern "C" double* TINKER_MOD(rgddyn, lm);
extern "C" double* TINKER_MOD(rgddyn, vc);
extern "C" double* TINKER_MOD(rgddyn, wc);
extern "C" int* TINKER_MOD(rgddyn, linear);

double*& xcmo = TINKER_MOD(rgddyn, xcmo);
double*& ycmo = TINKER_MOD(rgddyn, ycmo);
double*& zcmo = TINKER_MOD(rgddyn, zcmo);
double*& vcm = TINKER_MOD(rgddyn, vcm);
double*& wcm = TINKER_MOD(rgddyn, wcm);
double*& lm = TINKER_MOD(rgddyn, lm);
double*& vc = TINKER_MOD(rgddyn, vc);
double*& wc = TINKER_MOD(rgddyn, wc);
int*& linear = TINKER_MOD(rgddyn, linear);
#endif
} }
