#pragma once

#include "macro.hh"

namespace tinker { namespace freeze {
extern int& nrat;
extern int& nratx;
extern int*& iratx;
extern int*& kratx;
extern int*& irat;
extern double& rateps;
extern double*& krat;
extern int& use_rattle;
extern int*& ratimage;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(freeze, nrat);
extern "C" int TINKER_MOD(freeze, nratx);
extern "C" int* TINKER_MOD(freeze, iratx);
extern "C" int* TINKER_MOD(freeze, kratx);
extern "C" int* TINKER_MOD(freeze, irat);
extern "C" double TINKER_MOD(freeze, rateps);
extern "C" double* TINKER_MOD(freeze, krat);
extern "C" int TINKER_MOD(freeze, use_rattle);
extern "C" int* TINKER_MOD(freeze, ratimage);

int& nrat = TINKER_MOD(freeze, nrat);
int& nratx = TINKER_MOD(freeze, nratx);
int*& iratx = TINKER_MOD(freeze, iratx);
int*& kratx = TINKER_MOD(freeze, kratx);
int*& irat = TINKER_MOD(freeze, irat);
double& rateps = TINKER_MOD(freeze, rateps);
double*& krat = TINKER_MOD(freeze, krat);
int& use_rattle = TINKER_MOD(freeze, use_rattle);
int*& ratimage = TINKER_MOD(freeze, ratimage);
#endif
} }
