#pragma once

#include "macro.hh"

namespace tinker { namespace syntrn {
extern double& tpath;
extern double& ppath;
extern double*& xmin1;
extern double*& xmin2;
extern double*& xm;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(syntrn, tpath);
extern "C" double TINKER_MOD(syntrn, ppath);
extern "C" double* TINKER_MOD(syntrn, xmin1);
extern "C" double* TINKER_MOD(syntrn, xmin2);
extern "C" double* TINKER_MOD(syntrn, xm);

double& tpath = TINKER_MOD(syntrn, tpath);
double& ppath = TINKER_MOD(syntrn, ppath);
double*& xmin1 = TINKER_MOD(syntrn, xmin1);
double*& xmin2 = TINKER_MOD(syntrn, xmin2);
double*& xm = TINKER_MOD(syntrn, xm);
#endif
} }
