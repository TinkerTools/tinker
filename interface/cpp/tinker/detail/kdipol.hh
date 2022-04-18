#pragma once

#include "macro.hh"

namespace tinker { namespace kdipol {
extern int& maxnd;
extern int& maxnd5;
extern int& maxnd4;
extern int& maxnd3;
extern double*& dpl;
extern double*& dpl5;
extern double*& dpl4;
extern double*& dpl3;
extern double*& pos;
extern double*& pos5;
extern double*& pos4;
extern double*& pos3;
extern char (*&kd)[8];
extern char (*&kd5)[8];
extern char (*&kd4)[8];
extern char (*&kd3)[8];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(kdipol, maxnd);
extern "C" int TINKER_MOD(kdipol, maxnd5);
extern "C" int TINKER_MOD(kdipol, maxnd4);
extern "C" int TINKER_MOD(kdipol, maxnd3);
extern "C" double* TINKER_MOD(kdipol, dpl);
extern "C" double* TINKER_MOD(kdipol, dpl5);
extern "C" double* TINKER_MOD(kdipol, dpl4);
extern "C" double* TINKER_MOD(kdipol, dpl3);
extern "C" double* TINKER_MOD(kdipol, pos);
extern "C" double* TINKER_MOD(kdipol, pos5);
extern "C" double* TINKER_MOD(kdipol, pos4);
extern "C" double* TINKER_MOD(kdipol, pos3);
extern "C" char (*TINKER_MOD(kdipol, kd))[8];
extern "C" char (*TINKER_MOD(kdipol, kd5))[8];
extern "C" char (*TINKER_MOD(kdipol, kd4))[8];
extern "C" char (*TINKER_MOD(kdipol, kd3))[8];

int& maxnd = TINKER_MOD(kdipol, maxnd);
int& maxnd5 = TINKER_MOD(kdipol, maxnd5);
int& maxnd4 = TINKER_MOD(kdipol, maxnd4);
int& maxnd3 = TINKER_MOD(kdipol, maxnd3);
double*& dpl = TINKER_MOD(kdipol, dpl);
double*& dpl5 = TINKER_MOD(kdipol, dpl5);
double*& dpl4 = TINKER_MOD(kdipol, dpl4);
double*& dpl3 = TINKER_MOD(kdipol, dpl3);
double*& pos = TINKER_MOD(kdipol, pos);
double*& pos5 = TINKER_MOD(kdipol, pos5);
double*& pos4 = TINKER_MOD(kdipol, pos4);
double*& pos3 = TINKER_MOD(kdipol, pos3);
char (*&kd)[8] = TINKER_MOD(kdipol, kd);
char (*&kd5)[8] = TINKER_MOD(kdipol, kd5);
char (*&kd4)[8] = TINKER_MOD(kdipol, kd4);
char (*&kd3)[8] = TINKER_MOD(kdipol, kd3);
#endif
} }
