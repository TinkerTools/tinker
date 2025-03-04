#pragma once

#include "macro.hh"

namespace tinker { namespace tritor {
extern int& ntritor;
extern int*& itritor;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(tritor, ntritor);
extern "C" int* TINKER_MOD(tritor, itritor);

int& ntritor = TINKER_MOD(tritor, ntritor);
int*& itritor = TINKER_MOD(tritor, itritor);
#endif
} }
