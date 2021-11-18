#pragma once

#include "macro.hh"

namespace tinker { namespace limits {
extern double& vdwcut;
extern double& repcut;
extern double& dispcut;
extern double& chgcut;
extern double& dplcut;
extern double& mpolecut;
extern double& ctrncut;
extern double& vdwtaper;
extern double& reptaper;
extern double& disptaper;
extern double& chgtaper;
extern double& dpltaper;
extern double& mpoletaper;
extern double& ctrntaper;
extern double& ewaldcut;
extern double& dewaldcut;
extern double& usolvcut;
extern int& use_ewald;
extern int& use_dewald;
extern int& use_lights;
extern int& use_list;
extern int& use_vlist;
extern int& use_dlist;
extern int& use_clist;
extern int& use_mlist;
extern int& use_ulist;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(limits, vdwcut);
extern "C" double TINKER_MOD(limits, repcut);
extern "C" double TINKER_MOD(limits, dispcut);
extern "C" double TINKER_MOD(limits, chgcut);
extern "C" double TINKER_MOD(limits, dplcut);
extern "C" double TINKER_MOD(limits, mpolecut);
extern "C" double TINKER_MOD(limits, ctrncut);
extern "C" double TINKER_MOD(limits, vdwtaper);
extern "C" double TINKER_MOD(limits, reptaper);
extern "C" double TINKER_MOD(limits, disptaper);
extern "C" double TINKER_MOD(limits, chgtaper);
extern "C" double TINKER_MOD(limits, dpltaper);
extern "C" double TINKER_MOD(limits, mpoletaper);
extern "C" double TINKER_MOD(limits, ctrntaper);
extern "C" double TINKER_MOD(limits, ewaldcut);
extern "C" double TINKER_MOD(limits, dewaldcut);
extern "C" double TINKER_MOD(limits, usolvcut);
extern "C" int TINKER_MOD(limits, use_ewald);
extern "C" int TINKER_MOD(limits, use_dewald);
extern "C" int TINKER_MOD(limits, use_lights);
extern "C" int TINKER_MOD(limits, use_list);
extern "C" int TINKER_MOD(limits, use_vlist);
extern "C" int TINKER_MOD(limits, use_dlist);
extern "C" int TINKER_MOD(limits, use_clist);
extern "C" int TINKER_MOD(limits, use_mlist);
extern "C" int TINKER_MOD(limits, use_ulist);

double& vdwcut = TINKER_MOD(limits, vdwcut);
double& repcut = TINKER_MOD(limits, repcut);
double& dispcut = TINKER_MOD(limits, dispcut);
double& chgcut = TINKER_MOD(limits, chgcut);
double& dplcut = TINKER_MOD(limits, dplcut);
double& mpolecut = TINKER_MOD(limits, mpolecut);
double& ctrncut = TINKER_MOD(limits, ctrncut);
double& vdwtaper = TINKER_MOD(limits, vdwtaper);
double& reptaper = TINKER_MOD(limits, reptaper);
double& disptaper = TINKER_MOD(limits, disptaper);
double& chgtaper = TINKER_MOD(limits, chgtaper);
double& dpltaper = TINKER_MOD(limits, dpltaper);
double& mpoletaper = TINKER_MOD(limits, mpoletaper);
double& ctrntaper = TINKER_MOD(limits, ctrntaper);
double& ewaldcut = TINKER_MOD(limits, ewaldcut);
double& dewaldcut = TINKER_MOD(limits, dewaldcut);
double& usolvcut = TINKER_MOD(limits, usolvcut);
int& use_ewald = TINKER_MOD(limits, use_ewald);
int& use_dewald = TINKER_MOD(limits, use_dewald);
int& use_lights = TINKER_MOD(limits, use_lights);
int& use_list = TINKER_MOD(limits, use_list);
int& use_vlist = TINKER_MOD(limits, use_vlist);
int& use_dlist = TINKER_MOD(limits, use_dlist);
int& use_clist = TINKER_MOD(limits, use_clist);
int& use_mlist = TINKER_MOD(limits, use_mlist);
int& use_ulist = TINKER_MOD(limits, use_ulist);
#endif
} }
