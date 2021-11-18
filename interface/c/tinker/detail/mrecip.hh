#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern double TINKER_MOD(mrecip, vmxx);
extern double TINKER_MOD(mrecip, vmyy);
extern double TINKER_MOD(mrecip, vmzz);
extern double TINKER_MOD(mrecip, vmxy);
extern double TINKER_MOD(mrecip, vmxz);
extern double TINKER_MOD(mrecip, vmyz);
extern double* TINKER_MOD(mrecip, cmp);
extern double* TINKER_MOD(mrecip, fmp);
extern double* TINKER_MOD(mrecip, cphi);
extern double* TINKER_MOD(mrecip, fphi);
#ifdef __cplusplus
}
#endif
