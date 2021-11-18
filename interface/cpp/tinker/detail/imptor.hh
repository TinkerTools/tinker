#pragma once

#include "macro.hh"

namespace tinker { namespace imptor {
extern int& nitors;
extern int*& iitors;
extern double*& itors1;
extern double*& itors2;
extern double*& itors3;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(imptor, nitors);
extern "C" int* TINKER_MOD(imptor, iitors);
extern "C" double* TINKER_MOD(imptor, itors1);
extern "C" double* TINKER_MOD(imptor, itors2);
extern "C" double* TINKER_MOD(imptor, itors3);

int& nitors = TINKER_MOD(imptor, nitors);
int*& iitors = TINKER_MOD(imptor, iitors);
double*& itors1 = TINKER_MOD(imptor, itors1);
double*& itors2 = TINKER_MOD(imptor, itors2);
double*& itors3 = TINKER_MOD(imptor, itors3);
#endif
} }
