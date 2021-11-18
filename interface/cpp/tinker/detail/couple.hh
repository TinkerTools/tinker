#pragma once

#include "macro.hh"
#include "sizes.hh"

namespace tinker { namespace couple {
using namespace sizes;

extern int (&n12)[maxatm];
extern int*& n13;
extern int*& n14;
extern int*& n15;
extern int (&i12)[maxatm][maxval];
extern int*& i13;
extern int*& i14;
extern int*& i15;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(couple, n12)[maxatm];
extern "C" int* TINKER_MOD(couple, n13);
extern "C" int* TINKER_MOD(couple, n14);
extern "C" int* TINKER_MOD(couple, n15);
extern "C" int TINKER_MOD(couple, i12)[maxatm][maxval];
extern "C" int* TINKER_MOD(couple, i13);
extern "C" int* TINKER_MOD(couple, i14);
extern "C" int* TINKER_MOD(couple, i15);

int (&n12)[maxatm] = TINKER_MOD(couple, n12);
int*& n13 = TINKER_MOD(couple, n13);
int*& n14 = TINKER_MOD(couple, n14);
int*& n15 = TINKER_MOD(couple, n15);
int (&i12)[maxatm][maxval] = TINKER_MOD(couple, i12);
int*& i13 = TINKER_MOD(couple, i13);
int*& i14 = TINKER_MOD(couple, i14);
int*& i15 = TINKER_MOD(couple, i15);
#endif
} }
