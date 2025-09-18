#pragma once

#include "macro.hh"

namespace tinker { namespace mdstuf {
extern int& nfree;
extern int& irest;
extern int& bmnmix;
extern int& nrespa;
extern int& dorest;
extern char (&integrate)[11];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(mdstuf, nfree);
extern "C" int TINKER_MOD(mdstuf, irest);
extern "C" int TINKER_MOD(mdstuf, bmnmix);
extern "C" int TINKER_MOD(mdstuf, nrespa);
extern "C" int TINKER_MOD(mdstuf, dorest);
extern "C" char TINKER_MOD(mdstuf, integrate)[11];

int& nfree = TINKER_MOD(mdstuf, nfree);
int& irest = TINKER_MOD(mdstuf, irest);
int& bmnmix = TINKER_MOD(mdstuf, bmnmix);
int& nrespa = TINKER_MOD(mdstuf, nrespa);
int& dorest = TINKER_MOD(mdstuf, dorest);
char (&integrate)[11] = TINKER_MOD(mdstuf, integrate);
#endif
} }
