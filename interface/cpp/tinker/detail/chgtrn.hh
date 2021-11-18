#pragma once

#include "macro.hh"

namespace tinker { namespace chgtrn {
extern int& nct;
extern double*& chgct;
extern double*& dmpct;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(chgtrn, nct);
extern "C" double* TINKER_MOD(chgtrn, chgct);
extern "C" double* TINKER_MOD(chgtrn, dmpct);

int& nct = TINKER_MOD(chgtrn, nct);
double*& chgct = TINKER_MOD(chgtrn, chgct);
double*& dmpct = TINKER_MOD(chgtrn, dmpct);
#endif
} }
