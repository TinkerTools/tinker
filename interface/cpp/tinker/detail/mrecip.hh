#pragma once

#include "macro.hh"

namespace tinker { namespace mrecip {
extern double& vmxx;
extern double& vmyy;
extern double& vmzz;
extern double& vmxy;
extern double& vmxz;
extern double& vmyz;
extern double*& cmp;
extern double*& fmp;
extern double*& cphi;
extern double*& fphi;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(mrecip, vmxx);
extern "C" double TINKER_MOD(mrecip, vmyy);
extern "C" double TINKER_MOD(mrecip, vmzz);
extern "C" double TINKER_MOD(mrecip, vmxy);
extern "C" double TINKER_MOD(mrecip, vmxz);
extern "C" double TINKER_MOD(mrecip, vmyz);
extern "C" double* TINKER_MOD(mrecip, cmp);
extern "C" double* TINKER_MOD(mrecip, fmp);
extern "C" double* TINKER_MOD(mrecip, cphi);
extern "C" double* TINKER_MOD(mrecip, fphi);

double& vmxx = TINKER_MOD(mrecip, vmxx);
double& vmyy = TINKER_MOD(mrecip, vmyy);
double& vmzz = TINKER_MOD(mrecip, vmzz);
double& vmxy = TINKER_MOD(mrecip, vmxy);
double& vmxz = TINKER_MOD(mrecip, vmxz);
double& vmyz = TINKER_MOD(mrecip, vmyz);
double*& cmp = TINKER_MOD(mrecip, cmp);
double*& fmp = TINKER_MOD(mrecip, fmp);
double*& cphi = TINKER_MOD(mrecip, cphi);
double*& fphi = TINKER_MOD(mrecip, fphi);
#endif
} }
