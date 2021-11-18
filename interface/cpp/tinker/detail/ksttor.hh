#pragma once

#include "macro.hh"

namespace tinker { namespace ksttor {
const int maxnbt = 500;
extern double (&btcon)[maxnbt][9];
extern char (&kbt)[maxnbt][16];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(ksttor, btcon)[maxnbt][9];
extern "C" char TINKER_MOD(ksttor, kbt)[maxnbt][16];

double (&btcon)[maxnbt][9] = TINKER_MOD(ksttor, btcon);
char (&kbt)[maxnbt][16] = TINKER_MOD(ksttor, kbt);
#endif
} }
