#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern double TINKER_MOD(limits, vdwcut);
extern double TINKER_MOD(limits, repcut);
extern double TINKER_MOD(limits, dispcut);
extern double TINKER_MOD(limits, chgcut);
extern double TINKER_MOD(limits, dplcut);
extern double TINKER_MOD(limits, mpolecut);
extern double TINKER_MOD(limits, ctrncut);
extern double TINKER_MOD(limits, vdwtaper);
extern double TINKER_MOD(limits, reptaper);
extern double TINKER_MOD(limits, disptaper);
extern double TINKER_MOD(limits, chgtaper);
extern double TINKER_MOD(limits, dpltaper);
extern double TINKER_MOD(limits, mpoletaper);
extern double TINKER_MOD(limits, ctrntaper);
extern double TINKER_MOD(limits, ewaldcut);
extern double TINKER_MOD(limits, dewaldcut);
extern double TINKER_MOD(limits, usolvcut);
extern int TINKER_MOD(limits, use_ewald);
extern int TINKER_MOD(limits, use_dewald);
extern int TINKER_MOD(limits, use_lights);
extern int TINKER_MOD(limits, use_list);
extern int TINKER_MOD(limits, use_vlist);
extern int TINKER_MOD(limits, use_dlist);
extern int TINKER_MOD(limits, use_clist);
extern int TINKER_MOD(limits, use_mlist);
extern int TINKER_MOD(limits, use_ulist);
#ifdef __cplusplus
}
#endif
