#pragma once

#include "macro.hh"

namespace tinker { namespace vdw {
extern int& nvdw;
extern int*& ivdw;
extern int*& jvdw;
extern int*& ired;
extern double*& kred;
extern double*& xred;
extern double*& yred;
extern double*& zred;
extern double*& radmin;
extern double*& epsilon;
extern double*& radmin4;
extern double*& epsilon4;
extern double*& radhbnd;
extern double*& epshbnd;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(vdw, nvdw);
extern "C" int* TINKER_MOD(vdw, ivdw);
extern "C" int* TINKER_MOD(vdw, jvdw);
extern "C" int* TINKER_MOD(vdw, ired);
extern "C" double* TINKER_MOD(vdw, kred);
extern "C" double* TINKER_MOD(vdw, xred);
extern "C" double* TINKER_MOD(vdw, yred);
extern "C" double* TINKER_MOD(vdw, zred);
extern "C" double* TINKER_MOD(vdw, radmin);
extern "C" double* TINKER_MOD(vdw, epsilon);
extern "C" double* TINKER_MOD(vdw, radmin4);
extern "C" double* TINKER_MOD(vdw, epsilon4);
extern "C" double* TINKER_MOD(vdw, radhbnd);
extern "C" double* TINKER_MOD(vdw, epshbnd);

int& nvdw = TINKER_MOD(vdw, nvdw);
int*& ivdw = TINKER_MOD(vdw, ivdw);
int*& jvdw = TINKER_MOD(vdw, jvdw);
int*& ired = TINKER_MOD(vdw, ired);
double*& kred = TINKER_MOD(vdw, kred);
double*& xred = TINKER_MOD(vdw, xred);
double*& yred = TINKER_MOD(vdw, yred);
double*& zred = TINKER_MOD(vdw, zred);
double*& radmin = TINKER_MOD(vdw, radmin);
double*& epsilon = TINKER_MOD(vdw, epsilon);
double*& radmin4 = TINKER_MOD(vdw, radmin4);
double*& epsilon4 = TINKER_MOD(vdw, epsilon4);
double*& radhbnd = TINKER_MOD(vdw, radhbnd);
double*& epshbnd = TINKER_MOD(vdw, epshbnd);
#endif
} }
