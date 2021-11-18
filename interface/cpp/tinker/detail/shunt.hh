#pragma once

#include "macro.hh"

namespace tinker { namespace shunt {
extern double& off;
extern double& off2;
extern double& cut;
extern double& cut2;
extern double& c0;
extern double& c1;
extern double& c2;
extern double& c3;
extern double& c4;
extern double& c5;
extern double& f0;
extern double& f1;
extern double& f2;
extern double& f3;
extern double& f4;
extern double& f5;
extern double& f6;
extern double& f7;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(shunt, off);
extern "C" double TINKER_MOD(shunt, off2);
extern "C" double TINKER_MOD(shunt, cut);
extern "C" double TINKER_MOD(shunt, cut2);
extern "C" double TINKER_MOD(shunt, c0);
extern "C" double TINKER_MOD(shunt, c1);
extern "C" double TINKER_MOD(shunt, c2);
extern "C" double TINKER_MOD(shunt, c3);
extern "C" double TINKER_MOD(shunt, c4);
extern "C" double TINKER_MOD(shunt, c5);
extern "C" double TINKER_MOD(shunt, f0);
extern "C" double TINKER_MOD(shunt, f1);
extern "C" double TINKER_MOD(shunt, f2);
extern "C" double TINKER_MOD(shunt, f3);
extern "C" double TINKER_MOD(shunt, f4);
extern "C" double TINKER_MOD(shunt, f5);
extern "C" double TINKER_MOD(shunt, f6);
extern "C" double TINKER_MOD(shunt, f7);

double& off = TINKER_MOD(shunt, off);
double& off2 = TINKER_MOD(shunt, off2);
double& cut = TINKER_MOD(shunt, cut);
double& cut2 = TINKER_MOD(shunt, cut2);
double& c0 = TINKER_MOD(shunt, c0);
double& c1 = TINKER_MOD(shunt, c1);
double& c2 = TINKER_MOD(shunt, c2);
double& c3 = TINKER_MOD(shunt, c3);
double& c4 = TINKER_MOD(shunt, c4);
double& c5 = TINKER_MOD(shunt, c5);
double& f0 = TINKER_MOD(shunt, f0);
double& f1 = TINKER_MOD(shunt, f1);
double& f2 = TINKER_MOD(shunt, f2);
double& f3 = TINKER_MOD(shunt, f3);
double& f4 = TINKER_MOD(shunt, f4);
double& f5 = TINKER_MOD(shunt, f5);
double& f6 = TINKER_MOD(shunt, f6);
double& f7 = TINKER_MOD(shunt, f7);
#endif
} }
