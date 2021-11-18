#pragma once

#include "macro.hh"

namespace tinker { namespace valfit {
extern int& fit_bond;
extern int& fit_angle;
extern int& fit_strbnd;
extern int& fit_urey;
extern int& fit_opbend;
extern int& fit_tors;
extern int& fit_struct;
extern int& fit_force;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(valfit, fit_bond);
extern "C" int TINKER_MOD(valfit, fit_angle);
extern "C" int TINKER_MOD(valfit, fit_strbnd);
extern "C" int TINKER_MOD(valfit, fit_urey);
extern "C" int TINKER_MOD(valfit, fit_opbend);
extern "C" int TINKER_MOD(valfit, fit_tors);
extern "C" int TINKER_MOD(valfit, fit_struct);
extern "C" int TINKER_MOD(valfit, fit_force);

int& fit_bond = TINKER_MOD(valfit, fit_bond);
int& fit_angle = TINKER_MOD(valfit, fit_angle);
int& fit_strbnd = TINKER_MOD(valfit, fit_strbnd);
int& fit_urey = TINKER_MOD(valfit, fit_urey);
int& fit_opbend = TINKER_MOD(valfit, fit_opbend);
int& fit_tors = TINKER_MOD(valfit, fit_tors);
int& fit_struct = TINKER_MOD(valfit, fit_struct);
int& fit_force = TINKER_MOD(valfit, fit_force);
#endif
} }
