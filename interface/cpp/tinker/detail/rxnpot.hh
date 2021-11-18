#pragma once

#include "macro.hh"

namespace tinker { namespace rxnpot {
extern int& rfterms;
extern double& rfsize;
extern double& rfbulkd;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(rxnpot, rfterms);
extern "C" double TINKER_MOD(rxnpot, rfsize);
extern "C" double TINKER_MOD(rxnpot, rfbulkd);

int& rfterms = TINKER_MOD(rxnpot, rfterms);
double& rfsize = TINKER_MOD(rxnpot, rfsize);
double& rfbulkd = TINKER_MOD(rxnpot, rfbulkd);
#endif
} }
