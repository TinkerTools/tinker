#pragma once

#include "macro.hh"

namespace tinker { namespace orbits {
extern double*& qorb;
extern double*& worb;
extern double*& emorb;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(orbits, qorb);
extern "C" double* TINKER_MOD(orbits, worb);
extern "C" double* TINKER_MOD(orbits, emorb);

double*& qorb = TINKER_MOD(orbits, qorb);
double*& worb = TINKER_MOD(orbits, worb);
double*& emorb = TINKER_MOD(orbits, emorb);
#endif
} }
