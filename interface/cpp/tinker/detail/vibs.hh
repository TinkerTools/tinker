#pragma once

#include "macro.hh"

namespace tinker { namespace vibs {
extern double*& rho;
extern double*& rhok;
extern double*& rwork;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(vibs, rho);
extern "C" double* TINKER_MOD(vibs, rhok);
extern "C" double* TINKER_MOD(vibs, rwork);

double*& rho = TINKER_MOD(vibs, rho);
double*& rhok = TINKER_MOD(vibs, rhok);
double*& rwork = TINKER_MOD(vibs, rwork);
#endif
} }
