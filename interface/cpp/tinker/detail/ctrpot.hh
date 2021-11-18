#pragma once

#include "macro.hh"

namespace tinker { namespace ctrpot {
extern char (&ctrntyp)[8];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" char TINKER_MOD(ctrpot, ctrntyp)[8];

char (&ctrntyp)[8] = TINKER_MOD(ctrpot, ctrntyp);
#endif
} }
