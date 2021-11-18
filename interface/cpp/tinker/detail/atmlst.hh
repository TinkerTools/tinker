#pragma once

#include "macro.hh"

namespace tinker { namespace atmlst {
extern int*& bndlist;
extern int*& anglist;
extern int*& balist;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int* TINKER_MOD(atmlst, bndlist);
extern "C" int* TINKER_MOD(atmlst, anglist);
extern "C" int* TINKER_MOD(atmlst, balist);

int*& bndlist = TINKER_MOD(atmlst, bndlist);
int*& anglist = TINKER_MOD(atmlst, anglist);
int*& balist = TINKER_MOD(atmlst, balist);
#endif
} }
