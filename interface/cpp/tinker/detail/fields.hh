#pragma once

#include "macro.hh"

namespace tinker { namespace fields {
extern int*& biotyp;
extern char (&forcefield)[20];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int* TINKER_MOD(fields, biotyp);
extern "C" char TINKER_MOD(fields, forcefield)[20];

int*& biotyp = TINKER_MOD(fields, biotyp);
char (&forcefield)[20] = TINKER_MOD(fields, forcefield);
#endif
} }
