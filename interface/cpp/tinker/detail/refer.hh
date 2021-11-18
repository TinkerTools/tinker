#pragma once

#include "macro.hh"
#include "sizes.hh"

namespace tinker { namespace refer {
using namespace sizes;

extern int (&nref)[maxref];
extern int (&refltitle)[maxref];
extern int (&refleng)[maxref];
extern int*& reftyp;
extern int*& n12ref;
extern int*& i12ref;
extern double (&xboxref)[maxref];
extern double (&yboxref)[maxref];
extern double (&zboxref)[maxref];
extern double (&alpharef)[maxref];
extern double (&betaref)[maxref];
extern double (&gammaref)[maxref];
extern double*& xref;
extern double*& yref;
extern double*& zref;
extern char (*&refnam)[3];
extern char (&reffile)[maxref][240];
extern char (&reftitle)[maxref][240];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(refer, nref)[maxref];
extern "C" int TINKER_MOD(refer, refltitle)[maxref];
extern "C" int TINKER_MOD(refer, refleng)[maxref];
extern "C" int* TINKER_MOD(refer, reftyp);
extern "C" int* TINKER_MOD(refer, n12ref);
extern "C" int* TINKER_MOD(refer, i12ref);
extern "C" double TINKER_MOD(refer, xboxref)[maxref];
extern "C" double TINKER_MOD(refer, yboxref)[maxref];
extern "C" double TINKER_MOD(refer, zboxref)[maxref];
extern "C" double TINKER_MOD(refer, alpharef)[maxref];
extern "C" double TINKER_MOD(refer, betaref)[maxref];
extern "C" double TINKER_MOD(refer, gammaref)[maxref];
extern "C" double* TINKER_MOD(refer, xref);
extern "C" double* TINKER_MOD(refer, yref);
extern "C" double* TINKER_MOD(refer, zref);
extern "C" char (*TINKER_MOD(refer, refnam))[3];
extern "C" char TINKER_MOD(refer, reffile)[maxref][240];
extern "C" char TINKER_MOD(refer, reftitle)[maxref][240];

int (&nref)[maxref] = TINKER_MOD(refer, nref);
int (&refltitle)[maxref] = TINKER_MOD(refer, refltitle);
int (&refleng)[maxref] = TINKER_MOD(refer, refleng);
int*& reftyp = TINKER_MOD(refer, reftyp);
int*& n12ref = TINKER_MOD(refer, n12ref);
int*& i12ref = TINKER_MOD(refer, i12ref);
double (&xboxref)[maxref] = TINKER_MOD(refer, xboxref);
double (&yboxref)[maxref] = TINKER_MOD(refer, yboxref);
double (&zboxref)[maxref] = TINKER_MOD(refer, zboxref);
double (&alpharef)[maxref] = TINKER_MOD(refer, alpharef);
double (&betaref)[maxref] = TINKER_MOD(refer, betaref);
double (&gammaref)[maxref] = TINKER_MOD(refer, gammaref);
double*& xref = TINKER_MOD(refer, xref);
double*& yref = TINKER_MOD(refer, yref);
double*& zref = TINKER_MOD(refer, zref);
char (*&refnam)[3] = TINKER_MOD(refer, refnam);
char (&reffile)[maxref][240] = TINKER_MOD(refer, reffile);
char (&reftitle)[maxref][240] = TINKER_MOD(refer, reftitle);
#endif
} }
