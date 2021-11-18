#pragma once

#include "macro.hh"

namespace tinker { namespace pme {
extern int& nfft1;
extern int& nfft2;
extern int& nfft3;
extern int& nefft1;
extern int& nefft2;
extern int& nefft3;
extern int& ndfft1;
extern int& ndfft2;
extern int& ndfft3;
extern int& bsorder;
extern int& bseorder;
extern int& bsporder;
extern int& bsdorder;
extern int*& igrid;
extern double*& bsmod1;
extern double*& bsmod2;
extern double*& bsmod3;
extern double*& bsbuild;
extern double*& thetai1;
extern double*& thetai2;
extern double*& thetai3;
extern double*& qgrid;
extern double*& qfac;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(pme, nfft1);
extern "C" int TINKER_MOD(pme, nfft2);
extern "C" int TINKER_MOD(pme, nfft3);
extern "C" int TINKER_MOD(pme, nefft1);
extern "C" int TINKER_MOD(pme, nefft2);
extern "C" int TINKER_MOD(pme, nefft3);
extern "C" int TINKER_MOD(pme, ndfft1);
extern "C" int TINKER_MOD(pme, ndfft2);
extern "C" int TINKER_MOD(pme, ndfft3);
extern "C" int TINKER_MOD(pme, bsorder);
extern "C" int TINKER_MOD(pme, bseorder);
extern "C" int TINKER_MOD(pme, bsporder);
extern "C" int TINKER_MOD(pme, bsdorder);
extern "C" int* TINKER_MOD(pme, igrid);
extern "C" double* TINKER_MOD(pme, bsmod1);
extern "C" double* TINKER_MOD(pme, bsmod2);
extern "C" double* TINKER_MOD(pme, bsmod3);
extern "C" double* TINKER_MOD(pme, bsbuild);
extern "C" double* TINKER_MOD(pme, thetai1);
extern "C" double* TINKER_MOD(pme, thetai2);
extern "C" double* TINKER_MOD(pme, thetai3);
extern "C" double* TINKER_MOD(pme, qgrid);
extern "C" double* TINKER_MOD(pme, qfac);

int& nfft1 = TINKER_MOD(pme, nfft1);
int& nfft2 = TINKER_MOD(pme, nfft2);
int& nfft3 = TINKER_MOD(pme, nfft3);
int& nefft1 = TINKER_MOD(pme, nefft1);
int& nefft2 = TINKER_MOD(pme, nefft2);
int& nefft3 = TINKER_MOD(pme, nefft3);
int& ndfft1 = TINKER_MOD(pme, ndfft1);
int& ndfft2 = TINKER_MOD(pme, ndfft2);
int& ndfft3 = TINKER_MOD(pme, ndfft3);
int& bsorder = TINKER_MOD(pme, bsorder);
int& bseorder = TINKER_MOD(pme, bseorder);
int& bsporder = TINKER_MOD(pme, bsporder);
int& bsdorder = TINKER_MOD(pme, bsdorder);
int*& igrid = TINKER_MOD(pme, igrid);
double*& bsmod1 = TINKER_MOD(pme, bsmod1);
double*& bsmod2 = TINKER_MOD(pme, bsmod2);
double*& bsmod3 = TINKER_MOD(pme, bsmod3);
double*& bsbuild = TINKER_MOD(pme, bsbuild);
double*& thetai1 = TINKER_MOD(pme, thetai1);
double*& thetai2 = TINKER_MOD(pme, thetai2);
double*& thetai3 = TINKER_MOD(pme, thetai3);
double*& qgrid = TINKER_MOD(pme, qgrid);
double*& qfac = TINKER_MOD(pme, qfac);
#endif
} }
