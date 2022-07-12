#pragma once

#include "macro.hh"

namespace tinker { namespace expol {
extern int& nexpol;
extern double*& kpep;
extern double*& prepep;
extern double*& dmppep;
extern double*& polscale;
extern double*& invpolscale;
extern int*& lpep;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(expol, nexpol);
extern "C" double* TINKER_MOD(expol, kpep);
extern "C" double* TINKER_MOD(expol, prepep);
extern "C" double* TINKER_MOD(expol, dmppep);
extern "C" double* TINKER_MOD(expol, polscale);
extern "C" double* TINKER_MOD(expol, invpolscale);
extern "C" int* TINKER_MOD(expol, lpep);

int& nexpol = TINKER_MOD(expol, nexpol);
double*& kpep = TINKER_MOD(expol, kpep);
double*& prepep = TINKER_MOD(expol, prepep);
double*& dmppep = TINKER_MOD(expol, dmppep);
double*& polscale = TINKER_MOD(expol, polscale);
double*& invpolscale = TINKER_MOD(expol, invpolscale);
int*& lpep = TINKER_MOD(expol, lpep);
#endif
} }
