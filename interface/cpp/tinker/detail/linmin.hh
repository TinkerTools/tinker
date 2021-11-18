#pragma once

#include "macro.hh"

namespace tinker { namespace linmin {
extern int& intmax;
extern double& stpmin;
extern double& stpmax;
extern double& cappa;
extern double& slpmax;
extern double& angmax;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(linmin, intmax);
extern "C" double TINKER_MOD(linmin, stpmin);
extern "C" double TINKER_MOD(linmin, stpmax);
extern "C" double TINKER_MOD(linmin, cappa);
extern "C" double TINKER_MOD(linmin, slpmax);
extern "C" double TINKER_MOD(linmin, angmax);

int& intmax = TINKER_MOD(linmin, intmax);
double& stpmin = TINKER_MOD(linmin, stpmin);
double& stpmax = TINKER_MOD(linmin, stpmax);
double& cappa = TINKER_MOD(linmin, cappa);
double& slpmax = TINKER_MOD(linmin, slpmax);
double& angmax = TINKER_MOD(linmin, angmax);
#endif
} }
