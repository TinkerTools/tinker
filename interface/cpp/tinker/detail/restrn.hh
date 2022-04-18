#pragma once

#include "macro.hh"

namespace tinker { namespace restrn {
extern int& maxfix;
extern int& npfix;
extern int& ndfix;
extern int& nafix;
extern int& ntfix;
extern int& ngfix;
extern int& nchir;
extern int*& ipfix;
extern int*& kpfix;
extern int*& idfix;
extern int*& iafix;
extern int*& itfix;
extern int*& igfix;
extern int*& ichir;
extern double& depth;
extern double& width;
extern double& rflat;
extern double& rwall;
extern double*& xpfix;
extern double*& ypfix;
extern double*& zpfix;
extern double*& pfix;
extern double*& dfix;
extern double*& afix;
extern double*& tfix;
extern double*& gfix;
extern double*& chir;
extern int& use_basin;
extern int& use_wall;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(restrn, maxfix);
extern "C" int TINKER_MOD(restrn, npfix);
extern "C" int TINKER_MOD(restrn, ndfix);
extern "C" int TINKER_MOD(restrn, nafix);
extern "C" int TINKER_MOD(restrn, ntfix);
extern "C" int TINKER_MOD(restrn, ngfix);
extern "C" int TINKER_MOD(restrn, nchir);
extern "C" int* TINKER_MOD(restrn, ipfix);
extern "C" int* TINKER_MOD(restrn, kpfix);
extern "C" int* TINKER_MOD(restrn, idfix);
extern "C" int* TINKER_MOD(restrn, iafix);
extern "C" int* TINKER_MOD(restrn, itfix);
extern "C" int* TINKER_MOD(restrn, igfix);
extern "C" int* TINKER_MOD(restrn, ichir);
extern "C" double TINKER_MOD(restrn, depth);
extern "C" double TINKER_MOD(restrn, width);
extern "C" double TINKER_MOD(restrn, rflat);
extern "C" double TINKER_MOD(restrn, rwall);
extern "C" double* TINKER_MOD(restrn, xpfix);
extern "C" double* TINKER_MOD(restrn, ypfix);
extern "C" double* TINKER_MOD(restrn, zpfix);
extern "C" double* TINKER_MOD(restrn, pfix);
extern "C" double* TINKER_MOD(restrn, dfix);
extern "C" double* TINKER_MOD(restrn, afix);
extern "C" double* TINKER_MOD(restrn, tfix);
extern "C" double* TINKER_MOD(restrn, gfix);
extern "C" double* TINKER_MOD(restrn, chir);
extern "C" int TINKER_MOD(restrn, use_basin);
extern "C" int TINKER_MOD(restrn, use_wall);

int& maxfix = TINKER_MOD(restrn, maxfix);
int& npfix = TINKER_MOD(restrn, npfix);
int& ndfix = TINKER_MOD(restrn, ndfix);
int& nafix = TINKER_MOD(restrn, nafix);
int& ntfix = TINKER_MOD(restrn, ntfix);
int& ngfix = TINKER_MOD(restrn, ngfix);
int& nchir = TINKER_MOD(restrn, nchir);
int*& ipfix = TINKER_MOD(restrn, ipfix);
int*& kpfix = TINKER_MOD(restrn, kpfix);
int*& idfix = TINKER_MOD(restrn, idfix);
int*& iafix = TINKER_MOD(restrn, iafix);
int*& itfix = TINKER_MOD(restrn, itfix);
int*& igfix = TINKER_MOD(restrn, igfix);
int*& ichir = TINKER_MOD(restrn, ichir);
double& depth = TINKER_MOD(restrn, depth);
double& width = TINKER_MOD(restrn, width);
double& rflat = TINKER_MOD(restrn, rflat);
double& rwall = TINKER_MOD(restrn, rwall);
double*& xpfix = TINKER_MOD(restrn, xpfix);
double*& ypfix = TINKER_MOD(restrn, ypfix);
double*& zpfix = TINKER_MOD(restrn, zpfix);
double*& pfix = TINKER_MOD(restrn, pfix);
double*& dfix = TINKER_MOD(restrn, dfix);
double*& afix = TINKER_MOD(restrn, afix);
double*& tfix = TINKER_MOD(restrn, tfix);
double*& gfix = TINKER_MOD(restrn, gfix);
double*& chir = TINKER_MOD(restrn, chir);
int& use_basin = TINKER_MOD(restrn, use_basin);
int& use_wall = TINKER_MOD(restrn, use_wall);
#endif
} }
