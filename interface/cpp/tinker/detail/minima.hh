#pragma once

#include "macro.hh"

namespace tinker { namespace minima {
extern int& maxiter;
extern int& nextiter;
extern double& fctmin;
extern double& hguess;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(minima, maxiter);
extern "C" int TINKER_MOD(minima, nextiter);
extern "C" double TINKER_MOD(minima, fctmin);
extern "C" double TINKER_MOD(minima, hguess);

int& maxiter = TINKER_MOD(minima, maxiter);
int& nextiter = TINKER_MOD(minima, nextiter);
double& fctmin = TINKER_MOD(minima, fctmin);
double& hguess = TINKER_MOD(minima, hguess);
#endif
} }
