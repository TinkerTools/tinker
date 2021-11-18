#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern double TINKER_MOD(disgeo, vdwmax);
extern double TINKER_MOD(disgeo, compact);
extern double TINKER_MOD(disgeo, pathmax);
extern double* TINKER_MOD(disgeo, dbnd);
extern double* TINKER_MOD(disgeo, georad);
extern int TINKER_MOD(disgeo, use_invert);
extern int TINKER_MOD(disgeo, use_anneal);
#ifdef __cplusplus
}
#endif
