#pragma once

#include "macro.hh"

namespace tinker { namespace kangs {
extern int& maxna;
extern int& maxna5;
extern int& maxna4;
extern int& maxna3;
extern int& maxnap;
extern int& maxnaf;
extern double*& acon;
extern double*& acon5;
extern double*& acon4;
extern double*& acon3;
extern double*& aconp;
extern double*& aconf;
extern double*& ang;
extern double*& ang5;
extern double*& ang4;
extern double*& ang3;
extern double*& angp;
extern double*& angf;
extern char (*&ka)[12];
extern char (*&ka5)[12];
extern char (*&ka4)[12];
extern char (*&ka3)[12];
extern char (*&kap)[12];
extern char (*&kaf)[12];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(kangs, maxna);
extern "C" int TINKER_MOD(kangs, maxna5);
extern "C" int TINKER_MOD(kangs, maxna4);
extern "C" int TINKER_MOD(kangs, maxna3);
extern "C" int TINKER_MOD(kangs, maxnap);
extern "C" int TINKER_MOD(kangs, maxnaf);
extern "C" double* TINKER_MOD(kangs, acon);
extern "C" double* TINKER_MOD(kangs, acon5);
extern "C" double* TINKER_MOD(kangs, acon4);
extern "C" double* TINKER_MOD(kangs, acon3);
extern "C" double* TINKER_MOD(kangs, aconp);
extern "C" double* TINKER_MOD(kangs, aconf);
extern "C" double* TINKER_MOD(kangs, ang);
extern "C" double* TINKER_MOD(kangs, ang5);
extern "C" double* TINKER_MOD(kangs, ang4);
extern "C" double* TINKER_MOD(kangs, ang3);
extern "C" double* TINKER_MOD(kangs, angp);
extern "C" double* TINKER_MOD(kangs, angf);
extern "C" char (*TINKER_MOD(kangs, ka))[12];
extern "C" char (*TINKER_MOD(kangs, ka5))[12];
extern "C" char (*TINKER_MOD(kangs, ka4))[12];
extern "C" char (*TINKER_MOD(kangs, ka3))[12];
extern "C" char (*TINKER_MOD(kangs, kap))[12];
extern "C" char (*TINKER_MOD(kangs, kaf))[12];

int& maxna = TINKER_MOD(kangs, maxna);
int& maxna5 = TINKER_MOD(kangs, maxna5);
int& maxna4 = TINKER_MOD(kangs, maxna4);
int& maxna3 = TINKER_MOD(kangs, maxna3);
int& maxnap = TINKER_MOD(kangs, maxnap);
int& maxnaf = TINKER_MOD(kangs, maxnaf);
double*& acon = TINKER_MOD(kangs, acon);
double*& acon5 = TINKER_MOD(kangs, acon5);
double*& acon4 = TINKER_MOD(kangs, acon4);
double*& acon3 = TINKER_MOD(kangs, acon3);
double*& aconp = TINKER_MOD(kangs, aconp);
double*& aconf = TINKER_MOD(kangs, aconf);
double*& ang = TINKER_MOD(kangs, ang);
double*& ang5 = TINKER_MOD(kangs, ang5);
double*& ang4 = TINKER_MOD(kangs, ang4);
double*& ang3 = TINKER_MOD(kangs, ang3);
double*& angp = TINKER_MOD(kangs, angp);
double*& angf = TINKER_MOD(kangs, angf);
char (*&ka)[12] = TINKER_MOD(kangs, ka);
char (*&ka5)[12] = TINKER_MOD(kangs, ka5);
char (*&ka4)[12] = TINKER_MOD(kangs, ka4);
char (*&ka3)[12] = TINKER_MOD(kangs, ka3);
char (*&kap)[12] = TINKER_MOD(kangs, kap);
char (*&kaf)[12] = TINKER_MOD(kangs, kaf);
#endif
} }
