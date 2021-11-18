#pragma once

#include "macro.hh"

namespace tinker { namespace ktorsn {
const int maxnt = 2000;
const int maxnt5 = 500;
const int maxnt4 = 500;
extern double (&t1)[maxnt][2];
extern double (&t2)[maxnt][2];
extern double (&t3)[maxnt][2];
extern double (&t4)[maxnt][2];
extern double (&t5)[maxnt][2];
extern double (&t6)[maxnt][2];
extern double (&t15)[maxnt5][2];
extern double (&t25)[maxnt5][2];
extern double (&t35)[maxnt5][2];
extern double (&t45)[maxnt5][2];
extern double (&t55)[maxnt5][2];
extern double (&t65)[maxnt5][2];
extern double (&t14)[maxnt4][2];
extern double (&t24)[maxnt4][2];
extern double (&t34)[maxnt4][2];
extern double (&t44)[maxnt4][2];
extern double (&t54)[maxnt4][2];
extern double (&t64)[maxnt4][2];
extern char (&kt)[maxnt][16];
extern char (&kt5)[maxnt5][16];
extern char (&kt4)[maxnt4][16];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(ktorsn, t1)[maxnt][2];
extern "C" double TINKER_MOD(ktorsn, t2)[maxnt][2];
extern "C" double TINKER_MOD(ktorsn, t3)[maxnt][2];
extern "C" double TINKER_MOD(ktorsn, t4)[maxnt][2];
extern "C" double TINKER_MOD(ktorsn, t5)[maxnt][2];
extern "C" double TINKER_MOD(ktorsn, t6)[maxnt][2];
extern "C" double TINKER_MOD(ktorsn, t15)[maxnt5][2];
extern "C" double TINKER_MOD(ktorsn, t25)[maxnt5][2];
extern "C" double TINKER_MOD(ktorsn, t35)[maxnt5][2];
extern "C" double TINKER_MOD(ktorsn, t45)[maxnt5][2];
extern "C" double TINKER_MOD(ktorsn, t55)[maxnt5][2];
extern "C" double TINKER_MOD(ktorsn, t65)[maxnt5][2];
extern "C" double TINKER_MOD(ktorsn, t14)[maxnt4][2];
extern "C" double TINKER_MOD(ktorsn, t24)[maxnt4][2];
extern "C" double TINKER_MOD(ktorsn, t34)[maxnt4][2];
extern "C" double TINKER_MOD(ktorsn, t44)[maxnt4][2];
extern "C" double TINKER_MOD(ktorsn, t54)[maxnt4][2];
extern "C" double TINKER_MOD(ktorsn, t64)[maxnt4][2];
extern "C" char TINKER_MOD(ktorsn, kt)[maxnt][16];
extern "C" char TINKER_MOD(ktorsn, kt5)[maxnt5][16];
extern "C" char TINKER_MOD(ktorsn, kt4)[maxnt4][16];

double (&t1)[maxnt][2] = TINKER_MOD(ktorsn, t1);
double (&t2)[maxnt][2] = TINKER_MOD(ktorsn, t2);
double (&t3)[maxnt][2] = TINKER_MOD(ktorsn, t3);
double (&t4)[maxnt][2] = TINKER_MOD(ktorsn, t4);
double (&t5)[maxnt][2] = TINKER_MOD(ktorsn, t5);
double (&t6)[maxnt][2] = TINKER_MOD(ktorsn, t6);
double (&t15)[maxnt5][2] = TINKER_MOD(ktorsn, t15);
double (&t25)[maxnt5][2] = TINKER_MOD(ktorsn, t25);
double (&t35)[maxnt5][2] = TINKER_MOD(ktorsn, t35);
double (&t45)[maxnt5][2] = TINKER_MOD(ktorsn, t45);
double (&t55)[maxnt5][2] = TINKER_MOD(ktorsn, t55);
double (&t65)[maxnt5][2] = TINKER_MOD(ktorsn, t65);
double (&t14)[maxnt4][2] = TINKER_MOD(ktorsn, t14);
double (&t24)[maxnt4][2] = TINKER_MOD(ktorsn, t24);
double (&t34)[maxnt4][2] = TINKER_MOD(ktorsn, t34);
double (&t44)[maxnt4][2] = TINKER_MOD(ktorsn, t44);
double (&t54)[maxnt4][2] = TINKER_MOD(ktorsn, t54);
double (&t64)[maxnt4][2] = TINKER_MOD(ktorsn, t64);
char (&kt)[maxnt][16] = TINKER_MOD(ktorsn, kt);
char (&kt5)[maxnt5][16] = TINKER_MOD(ktorsn, kt5);
char (&kt4)[maxnt4][16] = TINKER_MOD(ktorsn, kt4);
#endif
} }
