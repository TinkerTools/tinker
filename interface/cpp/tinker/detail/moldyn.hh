#pragma once

#include "macro.hh"

namespace tinker { namespace moldyn {
extern double*& v;
extern double*& a;
extern double*& aalt;
extern double*& aslow;
extern double*& afast;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(moldyn, v);
extern "C" double* TINKER_MOD(moldyn, a);
extern "C" double* TINKER_MOD(moldyn, aalt);
extern "C" double* TINKER_MOD(moldyn, aslow);
extern "C" double* TINKER_MOD(moldyn, afast);

double*& v = TINKER_MOD(moldyn, v);
double*& a = TINKER_MOD(moldyn, a);
double*& aalt = TINKER_MOD(moldyn, aalt);
double*& aslow = TINKER_MOD(moldyn, aslow);
double*& afast = TINKER_MOD(moldyn, afast);
#endif
} }
