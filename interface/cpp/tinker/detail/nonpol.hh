#pragma once

#include "macro.hh"

namespace tinker { namespace nonpol {
const double epso = 0.1100e0;
const double epsh = 0.0135e0;
const double rmino = 1.7025e0;
const double rminh = 1.3275e0;
const double awater = 0.033428e0;
const double slevy = 1.0e0;
const double shctd = 0.75e0;
const double cavoff = 0.0e0;
const double dspoff = 1.056e0;
extern double& solvprs;
extern double& surften;
extern double& spcut;
extern double& spoff;
extern double& stcut;
extern double& stoff;
extern double*& radcav;
extern double*& raddsp;
extern double*& epsdsp;
extern double*& cdsp;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(nonpol, solvprs);
extern "C" double TINKER_MOD(nonpol, surften);
extern "C" double TINKER_MOD(nonpol, spcut);
extern "C" double TINKER_MOD(nonpol, spoff);
extern "C" double TINKER_MOD(nonpol, stcut);
extern "C" double TINKER_MOD(nonpol, stoff);
extern "C" double* TINKER_MOD(nonpol, radcav);
extern "C" double* TINKER_MOD(nonpol, raddsp);
extern "C" double* TINKER_MOD(nonpol, epsdsp);
extern "C" double* TINKER_MOD(nonpol, cdsp);

double& solvprs = TINKER_MOD(nonpol, solvprs);
double& surften = TINKER_MOD(nonpol, surften);
double& spcut = TINKER_MOD(nonpol, spcut);
double& spoff = TINKER_MOD(nonpol, spoff);
double& stcut = TINKER_MOD(nonpol, stcut);
double& stoff = TINKER_MOD(nonpol, stoff);
double*& radcav = TINKER_MOD(nonpol, radcav);
double*& raddsp = TINKER_MOD(nonpol, raddsp);
double*& epsdsp = TINKER_MOD(nonpol, epsdsp);
double*& cdsp = TINKER_MOD(nonpol, cdsp);
#endif
} }
