#pragma once

#include "macro.hh"

namespace tinker { namespace polopt {
const int maxopt = 6;
extern int& optorder;
extern int& optlevel;
extern double*& copt;
extern double*& copm;
extern double*& uopt;
extern double*& uoptp;
extern double*& uopts;
extern double*& uoptps;
extern double*& fopt;
extern double*& foptp;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(polopt, optorder);
extern "C" int TINKER_MOD(polopt, optlevel);
extern "C" double* TINKER_MOD(polopt, copt);
extern "C" double* TINKER_MOD(polopt, copm);
extern "C" double* TINKER_MOD(polopt, uopt);
extern "C" double* TINKER_MOD(polopt, uoptp);
extern "C" double* TINKER_MOD(polopt, uopts);
extern "C" double* TINKER_MOD(polopt, uoptps);
extern "C" double* TINKER_MOD(polopt, fopt);
extern "C" double* TINKER_MOD(polopt, foptp);

int& optorder = TINKER_MOD(polopt, optorder);
int& optlevel = TINKER_MOD(polopt, optlevel);
double*& copt = TINKER_MOD(polopt, copt);
double*& copm = TINKER_MOD(polopt, copm);
double*& uopt = TINKER_MOD(polopt, uopt);
double*& uoptp = TINKER_MOD(polopt, uoptp);
double*& uopts = TINKER_MOD(polopt, uopts);
double*& uoptps = TINKER_MOD(polopt, uoptps);
double*& fopt = TINKER_MOD(polopt, fopt);
double*& foptp = TINKER_MOD(polopt, foptp);
#endif
} }
