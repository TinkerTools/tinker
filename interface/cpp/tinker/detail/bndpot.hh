#pragma once

#include "macro.hh"

namespace tinker { namespace bndpot {
extern double& cbnd;
extern double& qbnd;
extern double& bndunit;
extern char (&bndtyp)[8];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(bndpot, cbnd);
extern "C" double TINKER_MOD(bndpot, qbnd);
extern "C" double TINKER_MOD(bndpot, bndunit);
extern "C" char TINKER_MOD(bndpot, bndtyp)[8];

double& cbnd = TINKER_MOD(bndpot, cbnd);
double& qbnd = TINKER_MOD(bndpot, qbnd);
double& bndunit = TINKER_MOD(bndpot, bndunit);
char (&bndtyp)[8] = TINKER_MOD(bndpot, bndtyp);
#endif
} }
