#pragma once

#include "macro.hh"
#include "sizes.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(nucleo, pucker)[TINKER_MOD__maxres];
extern double TINKER_MOD(nucleo, glyco)[TINKER_MOD__maxres];
extern double TINKER_MOD(nucleo, bkbone)[TINKER_MOD__maxres][6];
extern int TINKER_MOD(nucleo, dblhlx);
extern int TINKER_MOD(nucleo, deoxy)[TINKER_MOD__maxres];
extern char TINKER_MOD(nucleo, hlxform)[1];
#ifdef __cplusplus
}
#endif
