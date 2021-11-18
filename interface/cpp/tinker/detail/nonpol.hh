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
const double dispoff = 1.056e0;
extern double& solvprs;
extern double& surften;
extern double& spcut;
extern double& spoff;
extern double& stcut;
extern double& stoff;
extern double*& rcav;
extern double*& rdisp;
extern double*& cdisp;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(nonpol, solvprs);
extern "C" double TINKER_MOD(nonpol, surften);
extern "C" double TINKER_MOD(nonpol, spcut);
extern "C" double TINKER_MOD(nonpol, spoff);
extern "C" double TINKER_MOD(nonpol, stcut);
extern "C" double TINKER_MOD(nonpol, stoff);
extern "C" double* TINKER_MOD(nonpol, rcav);
extern "C" double* TINKER_MOD(nonpol, rdisp);
extern "C" double* TINKER_MOD(nonpol, cdisp);

double& solvprs = TINKER_MOD(nonpol, solvprs);
double& surften = TINKER_MOD(nonpol, surften);
double& spcut = TINKER_MOD(nonpol, spcut);
double& spoff = TINKER_MOD(nonpol, spoff);
double& stcut = TINKER_MOD(nonpol, stcut);
double& stoff = TINKER_MOD(nonpol, stoff);
double*& rcav = TINKER_MOD(nonpol, rcav);
double*& rdisp = TINKER_MOD(nonpol, rdisp);
double*& cdisp = TINKER_MOD(nonpol, cdisp);
#endif
} }
