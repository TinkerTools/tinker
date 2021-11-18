#pragma once

#include "macro.hh"

namespace tinker { namespace cflux {
extern int& nbflx;
extern int& naflx;
extern double*& bflx;
extern double*& aflx;
extern double*& abflx;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(cflux, nbflx);
extern "C" int TINKER_MOD(cflux, naflx);
extern "C" double* TINKER_MOD(cflux, bflx);
extern "C" double* TINKER_MOD(cflux, aflx);
extern "C" double* TINKER_MOD(cflux, abflx);

int& nbflx = TINKER_MOD(cflux, nbflx);
int& naflx = TINKER_MOD(cflux, naflx);
double*& bflx = TINKER_MOD(cflux, bflx);
double*& aflx = TINKER_MOD(cflux, aflx);
double*& abflx = TINKER_MOD(cflux, abflx);
#endif
} }
