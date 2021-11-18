#pragma once

#include "macro.hh"

namespace tinker { namespace tors {
extern int& ntors;
extern int*& itors;
extern double*& tors1;
extern double*& tors2;
extern double*& tors3;
extern double*& tors4;
extern double*& tors5;
extern double*& tors6;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(tors, ntors);
extern "C" int* TINKER_MOD(tors, itors);
extern "C" double* TINKER_MOD(tors, tors1);
extern "C" double* TINKER_MOD(tors, tors2);
extern "C" double* TINKER_MOD(tors, tors3);
extern "C" double* TINKER_MOD(tors, tors4);
extern "C" double* TINKER_MOD(tors, tors5);
extern "C" double* TINKER_MOD(tors, tors6);

int& ntors = TINKER_MOD(tors, ntors);
int*& itors = TINKER_MOD(tors, itors);
double*& tors1 = TINKER_MOD(tors, tors1);
double*& tors2 = TINKER_MOD(tors, tors2);
double*& tors3 = TINKER_MOD(tors, tors3);
double*& tors4 = TINKER_MOD(tors, tors4);
double*& tors5 = TINKER_MOD(tors, tors5);
double*& tors6 = TINKER_MOD(tors, tors6);
#endif
} }
