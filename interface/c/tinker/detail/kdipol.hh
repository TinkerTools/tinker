#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxnd 1000
#define TINKER_MOD__maxnd5 500
#define TINKER_MOD__maxnd4 500
#define TINKER_MOD__maxnd3 500
extern double TINKER_MOD(kdipol, dpl)[TINKER_MOD__maxnd];
extern double TINKER_MOD(kdipol, dpl5)[TINKER_MOD__maxnd5];
extern double TINKER_MOD(kdipol, dpl4)[TINKER_MOD__maxnd4];
extern double TINKER_MOD(kdipol, dpl3)[TINKER_MOD__maxnd3];
extern double TINKER_MOD(kdipol, pos)[TINKER_MOD__maxnd];
extern double TINKER_MOD(kdipol, pos5)[TINKER_MOD__maxnd5];
extern double TINKER_MOD(kdipol, pos4)[TINKER_MOD__maxnd4];
extern double TINKER_MOD(kdipol, pos3)[TINKER_MOD__maxnd3];
extern char TINKER_MOD(kdipol, kd)[TINKER_MOD__maxnd][8];
extern char TINKER_MOD(kdipol, kd5)[TINKER_MOD__maxnd5][8];
extern char TINKER_MOD(kdipol, kd4)[TINKER_MOD__maxnd4][8];
extern char TINKER_MOD(kdipol, kd3)[TINKER_MOD__maxnd3][8];
#ifdef __cplusplus
}
#endif
