#pragma once

#include "macro.hh"

namespace tinker { namespace energi {
extern double& esum;
extern double& eb;
extern double& ea;
extern double& eba;
extern double& eub;
extern double& eaa;
extern double& eopb;
extern double& eopd;
extern double& eid;
extern double& eit;
extern double& et;
extern double& ept;
extern double& ebt;
extern double& eat;
extern double& ett;
extern double& ev;
extern double& er;
extern double& edsp;
extern double& ec;
extern double& ecd;
extern double& ed;
extern double& em;
extern double& ep;
extern double& ect;
extern double& erxf;
extern double& es;
extern double& elf;
extern double& eg;
extern double& ex;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(energi, esum);
extern "C" double TINKER_MOD(energi, eb);
extern "C" double TINKER_MOD(energi, ea);
extern "C" double TINKER_MOD(energi, eba);
extern "C" double TINKER_MOD(energi, eub);
extern "C" double TINKER_MOD(energi, eaa);
extern "C" double TINKER_MOD(energi, eopb);
extern "C" double TINKER_MOD(energi, eopd);
extern "C" double TINKER_MOD(energi, eid);
extern "C" double TINKER_MOD(energi, eit);
extern "C" double TINKER_MOD(energi, et);
extern "C" double TINKER_MOD(energi, ept);
extern "C" double TINKER_MOD(energi, ebt);
extern "C" double TINKER_MOD(energi, eat);
extern "C" double TINKER_MOD(energi, ett);
extern "C" double TINKER_MOD(energi, ev);
extern "C" double TINKER_MOD(energi, er);
extern "C" double TINKER_MOD(energi, edsp);
extern "C" double TINKER_MOD(energi, ec);
extern "C" double TINKER_MOD(energi, ecd);
extern "C" double TINKER_MOD(energi, ed);
extern "C" double TINKER_MOD(energi, em);
extern "C" double TINKER_MOD(energi, ep);
extern "C" double TINKER_MOD(energi, ect);
extern "C" double TINKER_MOD(energi, erxf);
extern "C" double TINKER_MOD(energi, es);
extern "C" double TINKER_MOD(energi, elf);
extern "C" double TINKER_MOD(energi, eg);
extern "C" double TINKER_MOD(energi, ex);

double& esum = TINKER_MOD(energi, esum);
double& eb = TINKER_MOD(energi, eb);
double& ea = TINKER_MOD(energi, ea);
double& eba = TINKER_MOD(energi, eba);
double& eub = TINKER_MOD(energi, eub);
double& eaa = TINKER_MOD(energi, eaa);
double& eopb = TINKER_MOD(energi, eopb);
double& eopd = TINKER_MOD(energi, eopd);
double& eid = TINKER_MOD(energi, eid);
double& eit = TINKER_MOD(energi, eit);
double& et = TINKER_MOD(energi, et);
double& ept = TINKER_MOD(energi, ept);
double& ebt = TINKER_MOD(energi, ebt);
double& eat = TINKER_MOD(energi, eat);
double& ett = TINKER_MOD(energi, ett);
double& ev = TINKER_MOD(energi, ev);
double& er = TINKER_MOD(energi, er);
double& edsp = TINKER_MOD(energi, edsp);
double& ec = TINKER_MOD(energi, ec);
double& ecd = TINKER_MOD(energi, ecd);
double& ed = TINKER_MOD(energi, ed);
double& em = TINKER_MOD(energi, em);
double& ep = TINKER_MOD(energi, ep);
double& ect = TINKER_MOD(energi, ect);
double& erxf = TINKER_MOD(energi, erxf);
double& es = TINKER_MOD(energi, es);
double& elf = TINKER_MOD(energi, elf);
double& eg = TINKER_MOD(energi, eg);
double& ex = TINKER_MOD(energi, ex);
#endif
} }
