#pragma once

#include "macro.hh"

namespace tinker { namespace mpole {
const int maxpole = 13;
extern int& npole;
extern int*& ipole;
extern int*& polsiz;
extern int*& pollist;
extern int*& zaxis;
extern int*& xaxis;
extern int*& yaxis;
extern double*& pole;
extern double*& rpole;
extern double*& mono0;
extern char (*&polaxe)[8];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(mpole, npole);
extern "C" int* TINKER_MOD(mpole, ipole);
extern "C" int* TINKER_MOD(mpole, polsiz);
extern "C" int* TINKER_MOD(mpole, pollist);
extern "C" int* TINKER_MOD(mpole, zaxis);
extern "C" int* TINKER_MOD(mpole, xaxis);
extern "C" int* TINKER_MOD(mpole, yaxis);
extern "C" double* TINKER_MOD(mpole, pole);
extern "C" double* TINKER_MOD(mpole, rpole);
extern "C" double* TINKER_MOD(mpole, mono0);
extern "C" char (*TINKER_MOD(mpole, polaxe))[8];

int& npole = TINKER_MOD(mpole, npole);
int*& ipole = TINKER_MOD(mpole, ipole);
int*& polsiz = TINKER_MOD(mpole, polsiz);
int*& pollist = TINKER_MOD(mpole, pollist);
int*& zaxis = TINKER_MOD(mpole, zaxis);
int*& xaxis = TINKER_MOD(mpole, xaxis);
int*& yaxis = TINKER_MOD(mpole, yaxis);
double*& pole = TINKER_MOD(mpole, pole);
double*& rpole = TINKER_MOD(mpole, rpole);
double*& mono0 = TINKER_MOD(mpole, mono0);
char (*&polaxe)[8] = TINKER_MOD(mpole, polaxe);
#endif
} }
