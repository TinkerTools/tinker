#pragma once

#include "macro.hh"

namespace tinker { namespace domega {
extern double*& tesum;
extern double*& teb;
extern double*& tea;
extern double*& teba;
extern double*& teub;
extern double*& teaa;
extern double*& teopb;
extern double*& teopd;
extern double*& teid;
extern double*& teit;
extern double*& tet;
extern double*& tept;
extern double*& tebt;
extern double*& teat;
extern double*& tett;
extern double*& tev;
extern double*& ter;
extern double*& tedsp;
extern double*& tec;
extern double*& tecd;
extern double*& ted;
extern double*& tem;
extern double*& tep;
extern double*& tect;
extern double*& terxf;
extern double*& tes;
extern double*& telf;
extern double*& teg;
extern double*& tex;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(domega, tesum);
extern "C" double* TINKER_MOD(domega, teb);
extern "C" double* TINKER_MOD(domega, tea);
extern "C" double* TINKER_MOD(domega, teba);
extern "C" double* TINKER_MOD(domega, teub);
extern "C" double* TINKER_MOD(domega, teaa);
extern "C" double* TINKER_MOD(domega, teopb);
extern "C" double* TINKER_MOD(domega, teopd);
extern "C" double* TINKER_MOD(domega, teid);
extern "C" double* TINKER_MOD(domega, teit);
extern "C" double* TINKER_MOD(domega, tet);
extern "C" double* TINKER_MOD(domega, tept);
extern "C" double* TINKER_MOD(domega, tebt);
extern "C" double* TINKER_MOD(domega, teat);
extern "C" double* TINKER_MOD(domega, tett);
extern "C" double* TINKER_MOD(domega, tev);
extern "C" double* TINKER_MOD(domega, ter);
extern "C" double* TINKER_MOD(domega, tedsp);
extern "C" double* TINKER_MOD(domega, tec);
extern "C" double* TINKER_MOD(domega, tecd);
extern "C" double* TINKER_MOD(domega, ted);
extern "C" double* TINKER_MOD(domega, tem);
extern "C" double* TINKER_MOD(domega, tep);
extern "C" double* TINKER_MOD(domega, tect);
extern "C" double* TINKER_MOD(domega, terxf);
extern "C" double* TINKER_MOD(domega, tes);
extern "C" double* TINKER_MOD(domega, telf);
extern "C" double* TINKER_MOD(domega, teg);
extern "C" double* TINKER_MOD(domega, tex);

double*& tesum = TINKER_MOD(domega, tesum);
double*& teb = TINKER_MOD(domega, teb);
double*& tea = TINKER_MOD(domega, tea);
double*& teba = TINKER_MOD(domega, teba);
double*& teub = TINKER_MOD(domega, teub);
double*& teaa = TINKER_MOD(domega, teaa);
double*& teopb = TINKER_MOD(domega, teopb);
double*& teopd = TINKER_MOD(domega, teopd);
double*& teid = TINKER_MOD(domega, teid);
double*& teit = TINKER_MOD(domega, teit);
double*& tet = TINKER_MOD(domega, tet);
double*& tept = TINKER_MOD(domega, tept);
double*& tebt = TINKER_MOD(domega, tebt);
double*& teat = TINKER_MOD(domega, teat);
double*& tett = TINKER_MOD(domega, tett);
double*& tev = TINKER_MOD(domega, tev);
double*& ter = TINKER_MOD(domega, ter);
double*& tedsp = TINKER_MOD(domega, tedsp);
double*& tec = TINKER_MOD(domega, tec);
double*& tecd = TINKER_MOD(domega, tecd);
double*& ted = TINKER_MOD(domega, ted);
double*& tem = TINKER_MOD(domega, tem);
double*& tep = TINKER_MOD(domega, tep);
double*& tect = TINKER_MOD(domega, tect);
double*& terxf = TINKER_MOD(domega, terxf);
double*& tes = TINKER_MOD(domega, tes);
double*& telf = TINKER_MOD(domega, telf);
double*& teg = TINKER_MOD(domega, teg);
double*& tex = TINKER_MOD(domega, tex);
#endif
} }
