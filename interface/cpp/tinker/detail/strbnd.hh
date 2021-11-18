#pragma once

#include "macro.hh"

namespace tinker { namespace strbnd {
extern int& nstrbnd;
extern int*& isb;
extern double*& sbk;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(strbnd, nstrbnd);
extern "C" int* TINKER_MOD(strbnd, isb);
extern "C" double* TINKER_MOD(strbnd, sbk);

int& nstrbnd = TINKER_MOD(strbnd, nstrbnd);
int*& isb = TINKER_MOD(strbnd, isb);
double*& sbk = TINKER_MOD(strbnd, sbk);
#endif
} }
