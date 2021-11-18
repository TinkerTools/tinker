#pragma once

#include "macro.hh"

namespace tinker { namespace torpot {
extern double& idihunit;
extern double& itorunit;
extern double& torsunit;
extern double& ptorunit;
extern double& storunit;
extern double& atorunit;
extern double& ttorunit;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(torpot, idihunit);
extern "C" double TINKER_MOD(torpot, itorunit);
extern "C" double TINKER_MOD(torpot, torsunit);
extern "C" double TINKER_MOD(torpot, ptorunit);
extern "C" double TINKER_MOD(torpot, storunit);
extern "C" double TINKER_MOD(torpot, atorunit);
extern "C" double TINKER_MOD(torpot, ttorunit);

double& idihunit = TINKER_MOD(torpot, idihunit);
double& itorunit = TINKER_MOD(torpot, itorunit);
double& torsunit = TINKER_MOD(torpot, torsunit);
double& ptorunit = TINKER_MOD(torpot, ptorunit);
double& storunit = TINKER_MOD(torpot, storunit);
double& atorunit = TINKER_MOD(torpot, atorunit);
double& ttorunit = TINKER_MOD(torpot, ttorunit);
#endif
} }
