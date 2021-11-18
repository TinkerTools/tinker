#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(valfit, fit_bond);
extern int TINKER_MOD(valfit, fit_angle);
extern int TINKER_MOD(valfit, fit_strbnd);
extern int TINKER_MOD(valfit, fit_urey);
extern int TINKER_MOD(valfit, fit_opbend);
extern int TINKER_MOD(valfit, fit_tors);
extern int TINKER_MOD(valfit, fit_struct);
extern int TINKER_MOD(valfit, fit_force);
#ifdef __cplusplus
}
#endif
