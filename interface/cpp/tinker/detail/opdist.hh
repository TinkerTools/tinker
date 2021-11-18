#pragma once

#include "macro.hh"

namespace tinker { namespace opdist {
extern int& nopdist;
extern int*& iopd;
extern double*& opdk;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(opdist, nopdist);
extern "C" int* TINKER_MOD(opdist, iopd);
extern "C" double* TINKER_MOD(opdist, opdk);

int& nopdist = TINKER_MOD(opdist, nopdist);
int*& iopd = TINKER_MOD(opdist, iopd);
double*& opdk = TINKER_MOD(opdist, opdk);
#endif
} }
