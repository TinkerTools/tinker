#pragma once

#include "macro.hh"

namespace tinker { namespace repel {
extern int& nrep;
extern double*& sizpr;
extern double*& dmppr;
extern double*& elepr;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(repel, nrep);
extern "C" double* TINKER_MOD(repel, sizpr);
extern "C" double* TINKER_MOD(repel, dmppr);
extern "C" double* TINKER_MOD(repel, elepr);

int& nrep = TINKER_MOD(repel, nrep);
double*& sizpr = TINKER_MOD(repel, sizpr);
double*& dmppr = TINKER_MOD(repel, dmppr);
double*& elepr = TINKER_MOD(repel, elepr);
#endif
} }
