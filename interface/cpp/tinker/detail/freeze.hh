#pragma once

#include "macro.hh"

namespace tinker { namespace freeze {
extern int& nrat;
extern int& nratx;
extern int& nwat;
extern int& nwat4;
extern int*& iratx;
extern int*& kratx;
extern int*& irat;
extern int*& iwat;
extern int*& iwat4;
extern double& rateps;
extern double*& krat;
extern double*& kwat;
extern double*& kwat4;
extern int& use_freeze;
extern int*& frzimage;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(freeze, nrat);
extern "C" int TINKER_MOD(freeze, nratx);
extern "C" int TINKER_MOD(freeze, nwat);
extern "C" int TINKER_MOD(freeze, nwat4);
extern "C" int* TINKER_MOD(freeze, iratx);
extern "C" int* TINKER_MOD(freeze, kratx);
extern "C" int* TINKER_MOD(freeze, irat);
extern "C" int* TINKER_MOD(freeze, iwat);
extern "C" int* TINKER_MOD(freeze, iwat4);
extern "C" double TINKER_MOD(freeze, rateps);
extern "C" double* TINKER_MOD(freeze, krat);
extern "C" double* TINKER_MOD(freeze, kwat);
extern "C" double* TINKER_MOD(freeze, kwat4);
extern "C" int TINKER_MOD(freeze, use_freeze);
extern "C" int* TINKER_MOD(freeze, frzimage);

int& nrat = TINKER_MOD(freeze, nrat);
int& nratx = TINKER_MOD(freeze, nratx);
int& nwat = TINKER_MOD(freeze, nwat);
int& nwat4 = TINKER_MOD(freeze, nwat4);
int*& iratx = TINKER_MOD(freeze, iratx);
int*& kratx = TINKER_MOD(freeze, kratx);
int*& irat = TINKER_MOD(freeze, irat);
int*& iwat = TINKER_MOD(freeze, iwat);
int*& iwat4 = TINKER_MOD(freeze, iwat4);
double& rateps = TINKER_MOD(freeze, rateps);
double*& krat = TINKER_MOD(freeze, krat);
double*& kwat = TINKER_MOD(freeze, kwat);
double*& kwat4 = TINKER_MOD(freeze, kwat4);
int& use_freeze = TINKER_MOD(freeze, use_freeze);
int*& frzimage = TINKER_MOD(freeze, frzimage);
#endif
} }
