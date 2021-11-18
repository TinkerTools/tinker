#pragma once

#include "macro.hh"

namespace tinker { namespace hessn {
extern double*& hessx;
extern double*& hessy;
extern double*& hessz;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(hessn, hessx);
extern "C" double* TINKER_MOD(hessn, hessy);
extern "C" double* TINKER_MOD(hessn, hessz);

double*& hessx = TINKER_MOD(hessn, hessx);
double*& hessy = TINKER_MOD(hessn, hessy);
double*& hessz = TINKER_MOD(hessn, hessz);
#endif
} }
