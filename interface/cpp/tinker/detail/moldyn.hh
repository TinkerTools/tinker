#pragma once

#include "macro.hh"

namespace tinker { namespace moldyn {
extern double*& v;
extern double*& a;
extern double*& aalt;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(moldyn, v);
extern "C" double* TINKER_MOD(moldyn, a);
extern "C" double* TINKER_MOD(moldyn, aalt);

double*& v = TINKER_MOD(moldyn, v);
double*& a = TINKER_MOD(moldyn, a);
double*& aalt = TINKER_MOD(moldyn, aalt);
#endif
} }
