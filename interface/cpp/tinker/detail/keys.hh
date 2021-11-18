#pragma once

#include "macro.hh"

namespace tinker { namespace keys {
const int maxkey = 25000;
extern int& nkey;
extern char (&keyline)[maxkey][240];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(keys, nkey);
extern "C" char TINKER_MOD(keys, keyline)[maxkey][240];

int& nkey = TINKER_MOD(keys, nkey);
char (&keyline)[maxkey][240] = TINKER_MOD(keys, keyline);
#endif
} }
