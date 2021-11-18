#pragma once

#include "macro.hh"

namespace tinker { namespace argue {
const int maxarg = 20;
extern int& narg;
extern int (&listarg)[maxarg+1];
extern char (&arg)[maxarg+1][240];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(argue, narg);
extern "C" int TINKER_MOD(argue, listarg)[maxarg+1];
extern "C" char TINKER_MOD(argue, arg)[maxarg+1][240];

int& narg = TINKER_MOD(argue, narg);
int (&listarg)[maxarg+1] = TINKER_MOD(argue, listarg);
char (&arg)[maxarg+1][240] = TINKER_MOD(argue, arg);
#endif
} }
