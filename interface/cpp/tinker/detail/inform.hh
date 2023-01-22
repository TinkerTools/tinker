#pragma once

#include "macro.hh"

namespace tinker { namespace inform {
const int maxask = 5;
extern int& gpucard;
extern int& digits;
extern int& iprint;
extern int& iwrite;
extern int& isend;
extern int& verbose;
extern int& debug;
extern int& silent;
extern int& holdup;
extern int& abort;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(inform, gpucard);
extern "C" int TINKER_MOD(inform, digits);
extern "C" int TINKER_MOD(inform, iprint);
extern "C" int TINKER_MOD(inform, iwrite);
extern "C" int TINKER_MOD(inform, isend);
extern "C" int TINKER_MOD(inform, verbose);
extern "C" int TINKER_MOD(inform, debug);
extern "C" int TINKER_MOD(inform, silent);
extern "C" int TINKER_MOD(inform, holdup);
extern "C" int TINKER_MOD(inform, abort);

int& gpucard = TINKER_MOD(inform, gpucard);
int& digits = TINKER_MOD(inform, digits);
int& iprint = TINKER_MOD(inform, iprint);
int& iwrite = TINKER_MOD(inform, iwrite);
int& isend = TINKER_MOD(inform, isend);
int& verbose = TINKER_MOD(inform, verbose);
int& debug = TINKER_MOD(inform, debug);
int& silent = TINKER_MOD(inform, silent);
int& holdup = TINKER_MOD(inform, holdup);
int& abort = TINKER_MOD(inform, abort);
#endif
} }
