#pragma once

#include "macro.hh"

namespace tinker { namespace khbond {
extern int& maxnhb;
extern double*& radhb;
extern double*& epshb;
extern char (*&khb)[8];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(khbond, maxnhb);
extern "C" double* TINKER_MOD(khbond, radhb);
extern "C" double* TINKER_MOD(khbond, epshb);
extern "C" char (*TINKER_MOD(khbond, khb))[8];

int& maxnhb = TINKER_MOD(khbond, maxnhb);
double*& radhb = TINKER_MOD(khbond, radhb);
double*& epshb = TINKER_MOD(khbond, epshb);
char (*&khb)[8] = TINKER_MOD(khbond, khb);
#endif
} }
