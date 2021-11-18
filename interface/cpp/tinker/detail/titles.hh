#pragma once

#include "macro.hh"

namespace tinker { namespace titles {
extern int& ltitle;
extern char (&title)[240];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(titles, ltitle);
extern "C" char TINKER_MOD(titles, title)[240];

int& ltitle = TINKER_MOD(titles, ltitle);
char (&title)[240] = TINKER_MOD(titles, title);
#endif
} }
