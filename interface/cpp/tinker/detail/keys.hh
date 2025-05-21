#pragma once

#include "macro.hh"

namespace tinker { namespace keys {
extern int& nkey;
extern char (*&keyline)[240];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(keys, nkey);
extern "C" char (*TINKER_MOD(keys, keyline))[240];

int& nkey = TINKER_MOD(keys, nkey);
char (*&keyline)[240] = TINKER_MOD(keys, keyline);
#endif
} }
