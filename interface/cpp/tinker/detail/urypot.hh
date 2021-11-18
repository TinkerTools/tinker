#pragma once

#include "macro.hh"

namespace tinker { namespace urypot {
extern double& cury;
extern double& qury;
extern double& ureyunit;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(urypot, cury);
extern "C" double TINKER_MOD(urypot, qury);
extern "C" double TINKER_MOD(urypot, ureyunit);

double& cury = TINKER_MOD(urypot, cury);
double& qury = TINKER_MOD(urypot, qury);
double& ureyunit = TINKER_MOD(urypot, ureyunit);
#endif
} }
