#pragma once

#include "macro.hh"

namespace tinker { namespace piorbs {
extern int& norbit;
extern int& nconj;
extern int& reorbit;
extern int& nbpi;
extern int& ntpi;
extern int*& iorbit;
extern int*& iconj;
extern int*& kconj;
extern int*& piperp;
extern int*& ibpi;
extern int*& itpi;
extern double*& pbpl;
extern double*& pnpl;
extern int*& listpi;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(piorbs, norbit);
extern "C" int TINKER_MOD(piorbs, nconj);
extern "C" int TINKER_MOD(piorbs, reorbit);
extern "C" int TINKER_MOD(piorbs, nbpi);
extern "C" int TINKER_MOD(piorbs, ntpi);
extern "C" int* TINKER_MOD(piorbs, iorbit);
extern "C" int* TINKER_MOD(piorbs, iconj);
extern "C" int* TINKER_MOD(piorbs, kconj);
extern "C" int* TINKER_MOD(piorbs, piperp);
extern "C" int* TINKER_MOD(piorbs, ibpi);
extern "C" int* TINKER_MOD(piorbs, itpi);
extern "C" double* TINKER_MOD(piorbs, pbpl);
extern "C" double* TINKER_MOD(piorbs, pnpl);
extern "C" int* TINKER_MOD(piorbs, listpi);

int& norbit = TINKER_MOD(piorbs, norbit);
int& nconj = TINKER_MOD(piorbs, nconj);
int& reorbit = TINKER_MOD(piorbs, reorbit);
int& nbpi = TINKER_MOD(piorbs, nbpi);
int& ntpi = TINKER_MOD(piorbs, ntpi);
int*& iorbit = TINKER_MOD(piorbs, iorbit);
int*& iconj = TINKER_MOD(piorbs, iconj);
int*& kconj = TINKER_MOD(piorbs, kconj);
int*& piperp = TINKER_MOD(piorbs, piperp);
int*& ibpi = TINKER_MOD(piorbs, ibpi);
int*& itpi = TINKER_MOD(piorbs, itpi);
double*& pbpl = TINKER_MOD(piorbs, pbpl);
double*& pnpl = TINKER_MOD(piorbs, pnpl);
int*& listpi = TINKER_MOD(piorbs, listpi);
#endif
} }
