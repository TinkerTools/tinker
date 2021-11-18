#pragma once

#include "macro.hh"

namespace tinker { namespace angpot {
extern double& angunit;
extern double& stbnunit;
extern double& aaunit;
extern double& opbunit;
extern double& opdunit;
extern double& cang;
extern double& qang;
extern double& pang;
extern double& sang;
extern double& copb;
extern double& qopb;
extern double& popb;
extern double& sopb;
extern double& copd;
extern double& qopd;
extern double& popd;
extern double& sopd;
extern char (&opbtyp)[8];
extern char (*&angtyp)[8];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(angpot, angunit);
extern "C" double TINKER_MOD(angpot, stbnunit);
extern "C" double TINKER_MOD(angpot, aaunit);
extern "C" double TINKER_MOD(angpot, opbunit);
extern "C" double TINKER_MOD(angpot, opdunit);
extern "C" double TINKER_MOD(angpot, cang);
extern "C" double TINKER_MOD(angpot, qang);
extern "C" double TINKER_MOD(angpot, pang);
extern "C" double TINKER_MOD(angpot, sang);
extern "C" double TINKER_MOD(angpot, copb);
extern "C" double TINKER_MOD(angpot, qopb);
extern "C" double TINKER_MOD(angpot, popb);
extern "C" double TINKER_MOD(angpot, sopb);
extern "C" double TINKER_MOD(angpot, copd);
extern "C" double TINKER_MOD(angpot, qopd);
extern "C" double TINKER_MOD(angpot, popd);
extern "C" double TINKER_MOD(angpot, sopd);
extern "C" char TINKER_MOD(angpot, opbtyp)[8];
extern "C" char (*TINKER_MOD(angpot, angtyp))[8];

double& angunit = TINKER_MOD(angpot, angunit);
double& stbnunit = TINKER_MOD(angpot, stbnunit);
double& aaunit = TINKER_MOD(angpot, aaunit);
double& opbunit = TINKER_MOD(angpot, opbunit);
double& opdunit = TINKER_MOD(angpot, opdunit);
double& cang = TINKER_MOD(angpot, cang);
double& qang = TINKER_MOD(angpot, qang);
double& pang = TINKER_MOD(angpot, pang);
double& sang = TINKER_MOD(angpot, sang);
double& copb = TINKER_MOD(angpot, copb);
double& qopb = TINKER_MOD(angpot, qopb);
double& popb = TINKER_MOD(angpot, popb);
double& sopb = TINKER_MOD(angpot, sopb);
double& copd = TINKER_MOD(angpot, copd);
double& qopd = TINKER_MOD(angpot, qopd);
double& popd = TINKER_MOD(angpot, popd);
double& sopd = TINKER_MOD(angpot, sopd);
char (&opbtyp)[8] = TINKER_MOD(angpot, opbtyp);
char (*&angtyp)[8] = TINKER_MOD(angpot, angtyp);
#endif
} }
