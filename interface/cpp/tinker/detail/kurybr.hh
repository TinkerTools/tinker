#pragma once

#include "macro.hh"

namespace tinker { namespace kurybr {
extern int& maxnu;
extern double*& ucon;
extern double*& dst13;
extern char (*&ku)[12];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(kurybr, maxnu);
extern "C" double* TINKER_MOD(kurybr, ucon);
extern "C" double* TINKER_MOD(kurybr, dst13);
extern "C" char (*TINKER_MOD(kurybr, ku))[12];

int& maxnu = TINKER_MOD(kurybr, maxnu);
double*& ucon = TINKER_MOD(kurybr, ucon);
double*& dst13 = TINKER_MOD(kurybr, dst13);
char (*&ku)[12] = TINKER_MOD(kurybr, ku);
#endif
} }
