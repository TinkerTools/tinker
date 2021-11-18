#pragma once

#include "macro.hh"

namespace tinker { namespace ewald {
extern double& aewald;
extern double& aeewald;
extern double& apewald;
extern double& adewald;
extern char (&boundary)[7];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(ewald, aewald);
extern "C" double TINKER_MOD(ewald, aeewald);
extern "C" double TINKER_MOD(ewald, apewald);
extern "C" double TINKER_MOD(ewald, adewald);
extern "C" char TINKER_MOD(ewald, boundary)[7];

double& aewald = TINKER_MOD(ewald, aewald);
double& aeewald = TINKER_MOD(ewald, aeewald);
double& apewald = TINKER_MOD(ewald, apewald);
double& adewald = TINKER_MOD(ewald, adewald);
char (&boundary)[7] = TINKER_MOD(ewald, boundary);
#endif
} }
