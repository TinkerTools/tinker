#pragma once

#include "macro.hh"

namespace tinker { namespace inter {
extern double& einter;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(inter, einter);

double& einter = TINKER_MOD(inter, einter);
#endif
} }
