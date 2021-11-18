#pragma once

#include "macro.hh"

namespace tinker { namespace paths {
extern double& pnorm;
extern double (&acoeff)[7][7];
extern double*& pc0;
extern double*& pc1;
extern double*& pvect;
extern double*& pstep;
extern double*& pzet;
extern double*& gc;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(paths, pnorm);
extern "C" double TINKER_MOD(paths, acoeff)[7][7];
extern "C" double* TINKER_MOD(paths, pc0);
extern "C" double* TINKER_MOD(paths, pc1);
extern "C" double* TINKER_MOD(paths, pvect);
extern "C" double* TINKER_MOD(paths, pstep);
extern "C" double* TINKER_MOD(paths, pzet);
extern "C" double* TINKER_MOD(paths, gc);

double& pnorm = TINKER_MOD(paths, pnorm);
double (&acoeff)[7][7] = TINKER_MOD(paths, acoeff);
double*& pc0 = TINKER_MOD(paths, pc0);
double*& pc1 = TINKER_MOD(paths, pc1);
double*& pvect = TINKER_MOD(paths, pvect);
double*& pstep = TINKER_MOD(paths, pstep);
double*& pzet = TINKER_MOD(paths, pzet);
double*& gc = TINKER_MOD(paths, gc);
#endif
} }
