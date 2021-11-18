#pragma once

#include "macro.hh"

namespace tinker { namespace bitor_ {
extern int& nbitor;
extern int*& ibitor;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(bitor, nbitor);
extern "C" int* TINKER_MOD(bitor, ibitor);

int& nbitor = TINKER_MOD(bitor, nbitor);
int*& ibitor = TINKER_MOD(bitor, ibitor);
#endif
} }
