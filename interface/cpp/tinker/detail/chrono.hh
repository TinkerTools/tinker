#pragma once

#include "macro.hh"

namespace tinker { namespace chrono {
extern double& twall;
extern double& tcpu;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(chrono, twall);
extern "C" double TINKER_MOD(chrono, tcpu);

double& twall = TINKER_MOD(chrono, twall);
double& tcpu = TINKER_MOD(chrono, tcpu);
#endif
} }
