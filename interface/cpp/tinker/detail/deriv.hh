#pragma once

#include "macro.hh"

namespace tinker { namespace deriv {
extern double*& desum;
extern double*& deb;
extern double*& dea;
extern double*& deba;
extern double*& deub;
extern double*& deaa;
extern double*& deopb;
extern double*& deopd;
extern double*& deid;
extern double*& deit;
extern double*& det;
extern double*& dept;
extern double*& debt;
extern double*& deat;
extern double*& dett;
extern double*& dev;
extern double*& der;
extern double*& dedsp;
extern double*& dec;
extern double*& decd;
extern double*& ded;
extern double*& dem;
extern double*& dep;
extern double*& dect;
extern double*& derxf;
extern double*& des;
extern double*& delf;
extern double*& deg;
extern double*& dex;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(deriv, desum);
extern "C" double* TINKER_MOD(deriv, deb);
extern "C" double* TINKER_MOD(deriv, dea);
extern "C" double* TINKER_MOD(deriv, deba);
extern "C" double* TINKER_MOD(deriv, deub);
extern "C" double* TINKER_MOD(deriv, deaa);
extern "C" double* TINKER_MOD(deriv, deopb);
extern "C" double* TINKER_MOD(deriv, deopd);
extern "C" double* TINKER_MOD(deriv, deid);
extern "C" double* TINKER_MOD(deriv, deit);
extern "C" double* TINKER_MOD(deriv, det);
extern "C" double* TINKER_MOD(deriv, dept);
extern "C" double* TINKER_MOD(deriv, debt);
extern "C" double* TINKER_MOD(deriv, deat);
extern "C" double* TINKER_MOD(deriv, dett);
extern "C" double* TINKER_MOD(deriv, dev);
extern "C" double* TINKER_MOD(deriv, der);
extern "C" double* TINKER_MOD(deriv, dedsp);
extern "C" double* TINKER_MOD(deriv, dec);
extern "C" double* TINKER_MOD(deriv, decd);
extern "C" double* TINKER_MOD(deriv, ded);
extern "C" double* TINKER_MOD(deriv, dem);
extern "C" double* TINKER_MOD(deriv, dep);
extern "C" double* TINKER_MOD(deriv, dect);
extern "C" double* TINKER_MOD(deriv, derxf);
extern "C" double* TINKER_MOD(deriv, des);
extern "C" double* TINKER_MOD(deriv, delf);
extern "C" double* TINKER_MOD(deriv, deg);
extern "C" double* TINKER_MOD(deriv, dex);

double*& desum = TINKER_MOD(deriv, desum);
double*& deb = TINKER_MOD(deriv, deb);
double*& dea = TINKER_MOD(deriv, dea);
double*& deba = TINKER_MOD(deriv, deba);
double*& deub = TINKER_MOD(deriv, deub);
double*& deaa = TINKER_MOD(deriv, deaa);
double*& deopb = TINKER_MOD(deriv, deopb);
double*& deopd = TINKER_MOD(deriv, deopd);
double*& deid = TINKER_MOD(deriv, deid);
double*& deit = TINKER_MOD(deriv, deit);
double*& det = TINKER_MOD(deriv, det);
double*& dept = TINKER_MOD(deriv, dept);
double*& debt = TINKER_MOD(deriv, debt);
double*& deat = TINKER_MOD(deriv, deat);
double*& dett = TINKER_MOD(deriv, dett);
double*& dev = TINKER_MOD(deriv, dev);
double*& der = TINKER_MOD(deriv, der);
double*& dedsp = TINKER_MOD(deriv, dedsp);
double*& dec = TINKER_MOD(deriv, dec);
double*& decd = TINKER_MOD(deriv, decd);
double*& ded = TINKER_MOD(deriv, ded);
double*& dem = TINKER_MOD(deriv, dem);
double*& dep = TINKER_MOD(deriv, dep);
double*& dect = TINKER_MOD(deriv, dect);
double*& derxf = TINKER_MOD(deriv, derxf);
double*& des = TINKER_MOD(deriv, des);
double*& delf = TINKER_MOD(deriv, delf);
double*& deg = TINKER_MOD(deriv, deg);
double*& dex = TINKER_MOD(deriv, dex);
#endif
} }
