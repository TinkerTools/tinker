#pragma once

#include "macro.hh"

namespace tinker { namespace ksttor {
extern int& maxnbt;
extern double*& btcon;
extern char (*&kbt)[16];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(ksttor, maxnbt);
extern "C" double* TINKER_MOD(ksttor, btcon);
extern "C" char (*TINKER_MOD(ksttor, kbt))[16];

int& maxnbt = TINKER_MOD(ksttor, maxnbt);
double*& btcon = TINKER_MOD(ksttor, btcon);
char (*&kbt)[16] = TINKER_MOD(ksttor, kbt);
#endif
} }
