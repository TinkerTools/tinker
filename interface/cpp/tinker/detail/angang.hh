#pragma once

#include "macro.hh"

namespace tinker { namespace angang {
extern int& nangang;
extern int*& iaa;
extern double*& kaa;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(angang, nangang);
extern "C" int* TINKER_MOD(angang, iaa);
extern "C" double* TINKER_MOD(angang, kaa);

int& nangang = TINKER_MOD(angang, nangang);
int*& iaa = TINKER_MOD(angang, iaa);
double*& kaa = TINKER_MOD(angang, kaa);
#endif
} }
