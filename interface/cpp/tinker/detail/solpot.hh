#pragma once

#include "macro.hh"

namespace tinker { namespace solpot {
extern char (&solvtyp)[8];
extern char (&borntyp)[8];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" char TINKER_MOD(solpot, solvtyp)[8];
extern "C" char TINKER_MOD(solpot, borntyp)[8];

char (&solvtyp)[8] = TINKER_MOD(solpot, solvtyp);
char (&borntyp)[8] = TINKER_MOD(solpot, borntyp);
#endif
} }
