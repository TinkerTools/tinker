#pragma once

#include "macro.hh"

namespace tinker { namespace kdsp {
extern double*& dspsix;
extern double*& dspdmp;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(kdsp, dspsix);
extern "C" double* TINKER_MOD(kdsp, dspdmp);

double*& dspsix = TINKER_MOD(kdsp, dspsix);
double*& dspdmp = TINKER_MOD(kdsp, dspdmp);
#endif
} }
