#pragma once

#include "macro.hh"

namespace tinker { namespace dma {
extern double*& mp;
extern double*& dpx;
extern double*& dpy;
extern double*& dpz;
extern double*& q20;
extern double*& q21c;
extern double*& q21s;
extern double*& q22c;
extern double*& q22s;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(dma, mp);
extern "C" double* TINKER_MOD(dma, dpx);
extern "C" double* TINKER_MOD(dma, dpy);
extern "C" double* TINKER_MOD(dma, dpz);
extern "C" double* TINKER_MOD(dma, q20);
extern "C" double* TINKER_MOD(dma, q21c);
extern "C" double* TINKER_MOD(dma, q21s);
extern "C" double* TINKER_MOD(dma, q22c);
extern "C" double* TINKER_MOD(dma, q22s);

double*& mp = TINKER_MOD(dma, mp);
double*& dpx = TINKER_MOD(dma, dpx);
double*& dpy = TINKER_MOD(dma, dpy);
double*& dpz = TINKER_MOD(dma, dpz);
double*& q20 = TINKER_MOD(dma, q20);
double*& q21c = TINKER_MOD(dma, q21c);
double*& q21s = TINKER_MOD(dma, q21s);
double*& q22c = TINKER_MOD(dma, q22c);
double*& q22s = TINKER_MOD(dma, q22s);
#endif
} }
