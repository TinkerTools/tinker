#pragma once

#include "macro.hh"

namespace tinker { namespace iounit {
extern int& input;
extern int& iout;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(iounit, input);
extern "C" int TINKER_MOD(iounit, iout);

int& input = TINKER_MOD(iounit, input);
int& iout = TINKER_MOD(iounit, iout);
#endif
} }
