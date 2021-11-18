#pragma once

#include "macro.hh"

namespace tinker { namespace ring {
extern int& nring3;
extern int& nring4;
extern int& nring5;
extern int& nring6;
extern int& nring7;
extern int*& iring3;
extern int*& iring4;
extern int*& iring5;
extern int*& iring6;
extern int*& iring7;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(ring, nring3);
extern "C" int TINKER_MOD(ring, nring4);
extern "C" int TINKER_MOD(ring, nring5);
extern "C" int TINKER_MOD(ring, nring6);
extern "C" int TINKER_MOD(ring, nring7);
extern "C" int* TINKER_MOD(ring, iring3);
extern "C" int* TINKER_MOD(ring, iring4);
extern "C" int* TINKER_MOD(ring, iring5);
extern "C" int* TINKER_MOD(ring, iring6);
extern "C" int* TINKER_MOD(ring, iring7);

int& nring3 = TINKER_MOD(ring, nring3);
int& nring4 = TINKER_MOD(ring, nring4);
int& nring5 = TINKER_MOD(ring, nring5);
int& nring6 = TINKER_MOD(ring, nring6);
int& nring7 = TINKER_MOD(ring, nring7);
int*& iring3 = TINKER_MOD(ring, iring3);
int*& iring4 = TINKER_MOD(ring, iring4);
int*& iring5 = TINKER_MOD(ring, iring5);
int*& iring6 = TINKER_MOD(ring, iring6);
int*& iring7 = TINKER_MOD(ring, iring7);
#endif
} }
