#pragma once

#include "macro.hh"

namespace tinker { namespace ktorsn {
extern int& maxnt;
extern int& maxnt5;
extern int& maxnt4;
extern double*& t1;
extern double*& t2;
extern double*& t3;
extern double*& t4;
extern double*& t5;
extern double*& t6;
extern double*& t15;
extern double*& t25;
extern double*& t35;
extern double*& t45;
extern double*& t55;
extern double*& t65;
extern double*& t14;
extern double*& t24;
extern double*& t34;
extern double*& t44;
extern double*& t54;
extern double*& t64;
extern char (*&kt)[16];
extern char (*&kt5)[16];
extern char (*&kt4)[16];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(ktorsn, maxnt);
extern "C" int TINKER_MOD(ktorsn, maxnt5);
extern "C" int TINKER_MOD(ktorsn, maxnt4);
extern "C" double* TINKER_MOD(ktorsn, t1);
extern "C" double* TINKER_MOD(ktorsn, t2);
extern "C" double* TINKER_MOD(ktorsn, t3);
extern "C" double* TINKER_MOD(ktorsn, t4);
extern "C" double* TINKER_MOD(ktorsn, t5);
extern "C" double* TINKER_MOD(ktorsn, t6);
extern "C" double* TINKER_MOD(ktorsn, t15);
extern "C" double* TINKER_MOD(ktorsn, t25);
extern "C" double* TINKER_MOD(ktorsn, t35);
extern "C" double* TINKER_MOD(ktorsn, t45);
extern "C" double* TINKER_MOD(ktorsn, t55);
extern "C" double* TINKER_MOD(ktorsn, t65);
extern "C" double* TINKER_MOD(ktorsn, t14);
extern "C" double* TINKER_MOD(ktorsn, t24);
extern "C" double* TINKER_MOD(ktorsn, t34);
extern "C" double* TINKER_MOD(ktorsn, t44);
extern "C" double* TINKER_MOD(ktorsn, t54);
extern "C" double* TINKER_MOD(ktorsn, t64);
extern "C" char (*TINKER_MOD(ktorsn, kt))[16];
extern "C" char (*TINKER_MOD(ktorsn, kt5))[16];
extern "C" char (*TINKER_MOD(ktorsn, kt4))[16];

int& maxnt = TINKER_MOD(ktorsn, maxnt);
int& maxnt5 = TINKER_MOD(ktorsn, maxnt5);
int& maxnt4 = TINKER_MOD(ktorsn, maxnt4);
double*& t1 = TINKER_MOD(ktorsn, t1);
double*& t2 = TINKER_MOD(ktorsn, t2);
double*& t3 = TINKER_MOD(ktorsn, t3);
double*& t4 = TINKER_MOD(ktorsn, t4);
double*& t5 = TINKER_MOD(ktorsn, t5);
double*& t6 = TINKER_MOD(ktorsn, t6);
double*& t15 = TINKER_MOD(ktorsn, t15);
double*& t25 = TINKER_MOD(ktorsn, t25);
double*& t35 = TINKER_MOD(ktorsn, t35);
double*& t45 = TINKER_MOD(ktorsn, t45);
double*& t55 = TINKER_MOD(ktorsn, t55);
double*& t65 = TINKER_MOD(ktorsn, t65);
double*& t14 = TINKER_MOD(ktorsn, t14);
double*& t24 = TINKER_MOD(ktorsn, t24);
double*& t34 = TINKER_MOD(ktorsn, t34);
double*& t44 = TINKER_MOD(ktorsn, t44);
double*& t54 = TINKER_MOD(ktorsn, t54);
double*& t64 = TINKER_MOD(ktorsn, t64);
char (*&kt)[16] = TINKER_MOD(ktorsn, kt);
char (*&kt5)[16] = TINKER_MOD(ktorsn, kt5);
char (*&kt4)[16] = TINKER_MOD(ktorsn, kt4);
#endif
} }
