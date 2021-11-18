#pragma once

#include "macro.hh"

namespace tinker { namespace ptable {
const int maxele = 112;
extern double (&atmass)[maxele];
extern double (&vdwrad)[maxele];
extern double (&covrad)[maxele];
extern char (&elemnt)[maxele][3];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(ptable, atmass)[maxele];
extern "C" double TINKER_MOD(ptable, vdwrad)[maxele];
extern "C" double TINKER_MOD(ptable, covrad)[maxele];
extern "C" char TINKER_MOD(ptable, elemnt)[maxele][3];

double (&atmass)[maxele] = TINKER_MOD(ptable, atmass);
double (&vdwrad)[maxele] = TINKER_MOD(ptable, vdwrad);
double (&covrad)[maxele] = TINKER_MOD(ptable, covrad);
char (&elemnt)[maxele][3] = TINKER_MOD(ptable, elemnt);
#endif
} }
