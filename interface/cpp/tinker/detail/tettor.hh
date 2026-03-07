#pragma once

#include "macro.hh"

namespace tinker { namespace tettor {
extern int& ntettor;
extern int*& itettor;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(tettor, ntettor);
extern "C" int* TINKER_MOD(tettor, itettor);

int& ntettor = TINKER_MOD(tettor, ntettor);
int*& itettor = TINKER_MOD(tettor, itettor);
#endif
} }
