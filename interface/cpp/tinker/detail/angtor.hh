#pragma once

#include "macro.hh"

namespace tinker { namespace angtor {
extern int& nangtor;
extern int*& iat;
extern double*& kant;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(angtor, nangtor);
extern "C" int* TINKER_MOD(angtor, iat);
extern "C" double* TINKER_MOD(angtor, kant);

int& nangtor = TINKER_MOD(angtor, nangtor);
int*& iat = TINKER_MOD(angtor, iat);
double*& kant = TINKER_MOD(angtor, kant);
#endif
} }
