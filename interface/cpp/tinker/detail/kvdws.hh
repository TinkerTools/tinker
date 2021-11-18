#pragma once

#include "macro.hh"

namespace tinker { namespace kvdws {
extern double*& rad;
extern double*& eps;
extern double*& rad4;
extern double*& eps4;
extern double*& reduct;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(kvdws, rad);
extern "C" double* TINKER_MOD(kvdws, eps);
extern "C" double* TINKER_MOD(kvdws, rad4);
extern "C" double* TINKER_MOD(kvdws, eps4);
extern "C" double* TINKER_MOD(kvdws, reduct);

double*& rad = TINKER_MOD(kvdws, rad);
double*& eps = TINKER_MOD(kvdws, eps);
double*& rad4 = TINKER_MOD(kvdws, rad4);
double*& eps4 = TINKER_MOD(kvdws, eps4);
double*& reduct = TINKER_MOD(kvdws, reduct);
#endif
} }
