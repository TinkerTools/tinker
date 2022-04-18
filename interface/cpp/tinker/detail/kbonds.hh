#pragma once

#include "macro.hh"

namespace tinker { namespace kbonds {
extern int& maxnb;
extern int& maxnb5;
extern int& maxnb4;
extern int& maxnb3;
extern int& maxnel;
extern double*& bcon;
extern double*& bcon5;
extern double*& bcon4;
extern double*& bcon3;
extern double*& blen;
extern double*& blen5;
extern double*& blen4;
extern double*& blen3;
extern double*& dlen;
extern char (*&kb)[8];
extern char (*&kb5)[8];
extern char (*&kb4)[8];
extern char (*&kb3)[8];
extern char (*&kel)[12];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(kbonds, maxnb);
extern "C" int TINKER_MOD(kbonds, maxnb5);
extern "C" int TINKER_MOD(kbonds, maxnb4);
extern "C" int TINKER_MOD(kbonds, maxnb3);
extern "C" int TINKER_MOD(kbonds, maxnel);
extern "C" double* TINKER_MOD(kbonds, bcon);
extern "C" double* TINKER_MOD(kbonds, bcon5);
extern "C" double* TINKER_MOD(kbonds, bcon4);
extern "C" double* TINKER_MOD(kbonds, bcon3);
extern "C" double* TINKER_MOD(kbonds, blen);
extern "C" double* TINKER_MOD(kbonds, blen5);
extern "C" double* TINKER_MOD(kbonds, blen4);
extern "C" double* TINKER_MOD(kbonds, blen3);
extern "C" double* TINKER_MOD(kbonds, dlen);
extern "C" char (*TINKER_MOD(kbonds, kb))[8];
extern "C" char (*TINKER_MOD(kbonds, kb5))[8];
extern "C" char (*TINKER_MOD(kbonds, kb4))[8];
extern "C" char (*TINKER_MOD(kbonds, kb3))[8];
extern "C" char (*TINKER_MOD(kbonds, kel))[12];

int& maxnb = TINKER_MOD(kbonds, maxnb);
int& maxnb5 = TINKER_MOD(kbonds, maxnb5);
int& maxnb4 = TINKER_MOD(kbonds, maxnb4);
int& maxnb3 = TINKER_MOD(kbonds, maxnb3);
int& maxnel = TINKER_MOD(kbonds, maxnel);
double*& bcon = TINKER_MOD(kbonds, bcon);
double*& bcon5 = TINKER_MOD(kbonds, bcon5);
double*& bcon4 = TINKER_MOD(kbonds, bcon4);
double*& bcon3 = TINKER_MOD(kbonds, bcon3);
double*& blen = TINKER_MOD(kbonds, blen);
double*& blen5 = TINKER_MOD(kbonds, blen5);
double*& blen4 = TINKER_MOD(kbonds, blen4);
double*& blen3 = TINKER_MOD(kbonds, blen3);
double*& dlen = TINKER_MOD(kbonds, dlen);
char (*&kb)[8] = TINKER_MOD(kbonds, kb);
char (*&kb5)[8] = TINKER_MOD(kbonds, kb5);
char (*&kb4)[8] = TINKER_MOD(kbonds, kb4);
char (*&kb3)[8] = TINKER_MOD(kbonds, kb3);
char (*&kel)[12] = TINKER_MOD(kbonds, kel);
#endif
} }
