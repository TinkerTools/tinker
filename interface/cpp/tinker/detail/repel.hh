#pragma once

#include "macro.hh"

namespace tinker { namespace repel {
extern int& nrep;
extern int*& irep;
extern double*& sizpr;
extern double*& dmppr;
extern double*& elepr;
extern double*& repole;
extern double*& rrepole;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(repel, nrep);
extern "C" int* TINKER_MOD(repel, irep);
extern "C" double* TINKER_MOD(repel, sizpr);
extern "C" double* TINKER_MOD(repel, dmppr);
extern "C" double* TINKER_MOD(repel, elepr);
extern "C" double* TINKER_MOD(repel, repole);
extern "C" double* TINKER_MOD(repel, rrepole);

int& nrep = TINKER_MOD(repel, nrep);
int*& irep = TINKER_MOD(repel, irep);
double*& sizpr = TINKER_MOD(repel, sizpr);
double*& dmppr = TINKER_MOD(repel, dmppr);
double*& elepr = TINKER_MOD(repel, elepr);
double*& repole = TINKER_MOD(repel, repole);
double*& rrepole = TINKER_MOD(repel, rrepole);
#endif
} }
