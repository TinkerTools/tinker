#pragma once

#include "macro.hh"

namespace tinker { namespace ksolut {
extern double*& pbr;
extern double*& csr;
extern double*& gkr;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(ksolut, pbr);
extern "C" double* TINKER_MOD(ksolut, csr);
extern "C" double* TINKER_MOD(ksolut, gkr);

double*& pbr = TINKER_MOD(ksolut, pbr);
double*& csr = TINKER_MOD(ksolut, csr);
double*& gkr = TINKER_MOD(ksolut, gkr);
#endif
} }
