#pragma once

#include "macro.hh"

namespace tinker { namespace kdipol {
const int maxnd = 1000;
const int maxnd5 = 500;
const int maxnd4 = 500;
const int maxnd3 = 500;
extern double (&dpl)[maxnd];
extern double (&dpl5)[maxnd5];
extern double (&dpl4)[maxnd4];
extern double (&dpl3)[maxnd3];
extern double (&pos)[maxnd];
extern double (&pos5)[maxnd5];
extern double (&pos4)[maxnd4];
extern double (&pos3)[maxnd3];
extern char (&kd)[maxnd][8];
extern char (&kd5)[maxnd5][8];
extern char (&kd4)[maxnd4][8];
extern char (&kd3)[maxnd3][8];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(kdipol, dpl)[maxnd];
extern "C" double TINKER_MOD(kdipol, dpl5)[maxnd5];
extern "C" double TINKER_MOD(kdipol, dpl4)[maxnd4];
extern "C" double TINKER_MOD(kdipol, dpl3)[maxnd3];
extern "C" double TINKER_MOD(kdipol, pos)[maxnd];
extern "C" double TINKER_MOD(kdipol, pos5)[maxnd5];
extern "C" double TINKER_MOD(kdipol, pos4)[maxnd4];
extern "C" double TINKER_MOD(kdipol, pos3)[maxnd3];
extern "C" char TINKER_MOD(kdipol, kd)[maxnd][8];
extern "C" char TINKER_MOD(kdipol, kd5)[maxnd5][8];
extern "C" char TINKER_MOD(kdipol, kd4)[maxnd4][8];
extern "C" char TINKER_MOD(kdipol, kd3)[maxnd3][8];

double (&dpl)[maxnd] = TINKER_MOD(kdipol, dpl);
double (&dpl5)[maxnd5] = TINKER_MOD(kdipol, dpl5);
double (&dpl4)[maxnd4] = TINKER_MOD(kdipol, dpl4);
double (&dpl3)[maxnd3] = TINKER_MOD(kdipol, dpl3);
double (&pos)[maxnd] = TINKER_MOD(kdipol, pos);
double (&pos5)[maxnd5] = TINKER_MOD(kdipol, pos5);
double (&pos4)[maxnd4] = TINKER_MOD(kdipol, pos4);
double (&pos3)[maxnd3] = TINKER_MOD(kdipol, pos3);
char (&kd)[maxnd][8] = TINKER_MOD(kdipol, kd);
char (&kd5)[maxnd5][8] = TINKER_MOD(kdipol, kd5);
char (&kd4)[maxnd4][8] = TINKER_MOD(kdipol, kd4);
char (&kd3)[maxnd3][8] = TINKER_MOD(kdipol, kd3);
#endif
} }
