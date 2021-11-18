#pragma once

#include "macro.hh"

namespace tinker { namespace analyz {
extern double*& aesum;
extern double*& aeb;
extern double*& aea;
extern double*& aeba;
extern double*& aeub;
extern double*& aeaa;
extern double*& aeopb;
extern double*& aeopd;
extern double*& aeid;
extern double*& aeit;
extern double*& aet;
extern double*& aept;
extern double*& aebt;
extern double*& aeat;
extern double*& aett;
extern double*& aev;
extern double*& aer;
extern double*& aedsp;
extern double*& aec;
extern double*& aecd;
extern double*& aed;
extern double*& aem;
extern double*& aep;
extern double*& aect;
extern double*& aerxf;
extern double*& aes;
extern double*& aelf;
extern double*& aeg;
extern double*& aex;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(analyz, aesum);
extern "C" double* TINKER_MOD(analyz, aeb);
extern "C" double* TINKER_MOD(analyz, aea);
extern "C" double* TINKER_MOD(analyz, aeba);
extern "C" double* TINKER_MOD(analyz, aeub);
extern "C" double* TINKER_MOD(analyz, aeaa);
extern "C" double* TINKER_MOD(analyz, aeopb);
extern "C" double* TINKER_MOD(analyz, aeopd);
extern "C" double* TINKER_MOD(analyz, aeid);
extern "C" double* TINKER_MOD(analyz, aeit);
extern "C" double* TINKER_MOD(analyz, aet);
extern "C" double* TINKER_MOD(analyz, aept);
extern "C" double* TINKER_MOD(analyz, aebt);
extern "C" double* TINKER_MOD(analyz, aeat);
extern "C" double* TINKER_MOD(analyz, aett);
extern "C" double* TINKER_MOD(analyz, aev);
extern "C" double* TINKER_MOD(analyz, aer);
extern "C" double* TINKER_MOD(analyz, aedsp);
extern "C" double* TINKER_MOD(analyz, aec);
extern "C" double* TINKER_MOD(analyz, aecd);
extern "C" double* TINKER_MOD(analyz, aed);
extern "C" double* TINKER_MOD(analyz, aem);
extern "C" double* TINKER_MOD(analyz, aep);
extern "C" double* TINKER_MOD(analyz, aect);
extern "C" double* TINKER_MOD(analyz, aerxf);
extern "C" double* TINKER_MOD(analyz, aes);
extern "C" double* TINKER_MOD(analyz, aelf);
extern "C" double* TINKER_MOD(analyz, aeg);
extern "C" double* TINKER_MOD(analyz, aex);

double*& aesum = TINKER_MOD(analyz, aesum);
double*& aeb = TINKER_MOD(analyz, aeb);
double*& aea = TINKER_MOD(analyz, aea);
double*& aeba = TINKER_MOD(analyz, aeba);
double*& aeub = TINKER_MOD(analyz, aeub);
double*& aeaa = TINKER_MOD(analyz, aeaa);
double*& aeopb = TINKER_MOD(analyz, aeopb);
double*& aeopd = TINKER_MOD(analyz, aeopd);
double*& aeid = TINKER_MOD(analyz, aeid);
double*& aeit = TINKER_MOD(analyz, aeit);
double*& aet = TINKER_MOD(analyz, aet);
double*& aept = TINKER_MOD(analyz, aept);
double*& aebt = TINKER_MOD(analyz, aebt);
double*& aeat = TINKER_MOD(analyz, aeat);
double*& aett = TINKER_MOD(analyz, aett);
double*& aev = TINKER_MOD(analyz, aev);
double*& aer = TINKER_MOD(analyz, aer);
double*& aedsp = TINKER_MOD(analyz, aedsp);
double*& aec = TINKER_MOD(analyz, aec);
double*& aecd = TINKER_MOD(analyz, aecd);
double*& aed = TINKER_MOD(analyz, aed);
double*& aem = TINKER_MOD(analyz, aem);
double*& aep = TINKER_MOD(analyz, aep);
double*& aect = TINKER_MOD(analyz, aect);
double*& aerxf = TINKER_MOD(analyz, aerxf);
double*& aes = TINKER_MOD(analyz, aes);
double*& aelf = TINKER_MOD(analyz, aelf);
double*& aeg = TINKER_MOD(analyz, aeg);
double*& aex = TINKER_MOD(analyz, aex);
#endif
} }
