#pragma once

#include "macro.hh"

namespace tinker { namespace potent {
extern int& use_bond;
extern int& use_angle;
extern int& use_strbnd;
extern int& use_urey;
extern int& use_angang;
extern int& use_opbend;
extern int& use_opdist;
extern int& use_improp;
extern int& use_imptor;
extern int& use_tors;
extern int& use_pitors;
extern int& use_strtor;
extern int& use_angtor;
extern int& use_tortor;
extern int& use_vdw;
extern int& use_repuls;
extern int& use_disp;
extern int& use_charge;
extern int& use_chgdpl;
extern int& use_dipole;
extern int& use_mpole;
extern int& use_polar;
extern int& use_chgtrn;
extern int& use_chgflx;
extern int& use_rxnfld;
extern int& use_solv;
extern int& use_metal;
extern int& use_geom;
extern int& use_extra;
extern int& use_born;
extern int& use_orbit;
extern int& use_mutate;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(potent, use_bond);
extern "C" int TINKER_MOD(potent, use_angle);
extern "C" int TINKER_MOD(potent, use_strbnd);
extern "C" int TINKER_MOD(potent, use_urey);
extern "C" int TINKER_MOD(potent, use_angang);
extern "C" int TINKER_MOD(potent, use_opbend);
extern "C" int TINKER_MOD(potent, use_opdist);
extern "C" int TINKER_MOD(potent, use_improp);
extern "C" int TINKER_MOD(potent, use_imptor);
extern "C" int TINKER_MOD(potent, use_tors);
extern "C" int TINKER_MOD(potent, use_pitors);
extern "C" int TINKER_MOD(potent, use_strtor);
extern "C" int TINKER_MOD(potent, use_angtor);
extern "C" int TINKER_MOD(potent, use_tortor);
extern "C" int TINKER_MOD(potent, use_vdw);
extern "C" int TINKER_MOD(potent, use_repuls);
extern "C" int TINKER_MOD(potent, use_disp);
extern "C" int TINKER_MOD(potent, use_charge);
extern "C" int TINKER_MOD(potent, use_chgdpl);
extern "C" int TINKER_MOD(potent, use_dipole);
extern "C" int TINKER_MOD(potent, use_mpole);
extern "C" int TINKER_MOD(potent, use_polar);
extern "C" int TINKER_MOD(potent, use_chgtrn);
extern "C" int TINKER_MOD(potent, use_chgflx);
extern "C" int TINKER_MOD(potent, use_rxnfld);
extern "C" int TINKER_MOD(potent, use_solv);
extern "C" int TINKER_MOD(potent, use_metal);
extern "C" int TINKER_MOD(potent, use_geom);
extern "C" int TINKER_MOD(potent, use_extra);
extern "C" int TINKER_MOD(potent, use_born);
extern "C" int TINKER_MOD(potent, use_orbit);
extern "C" int TINKER_MOD(potent, use_mutate);

int& use_bond = TINKER_MOD(potent, use_bond);
int& use_angle = TINKER_MOD(potent, use_angle);
int& use_strbnd = TINKER_MOD(potent, use_strbnd);
int& use_urey = TINKER_MOD(potent, use_urey);
int& use_angang = TINKER_MOD(potent, use_angang);
int& use_opbend = TINKER_MOD(potent, use_opbend);
int& use_opdist = TINKER_MOD(potent, use_opdist);
int& use_improp = TINKER_MOD(potent, use_improp);
int& use_imptor = TINKER_MOD(potent, use_imptor);
int& use_tors = TINKER_MOD(potent, use_tors);
int& use_pitors = TINKER_MOD(potent, use_pitors);
int& use_strtor = TINKER_MOD(potent, use_strtor);
int& use_angtor = TINKER_MOD(potent, use_angtor);
int& use_tortor = TINKER_MOD(potent, use_tortor);
int& use_vdw = TINKER_MOD(potent, use_vdw);
int& use_repuls = TINKER_MOD(potent, use_repuls);
int& use_disp = TINKER_MOD(potent, use_disp);
int& use_charge = TINKER_MOD(potent, use_charge);
int& use_chgdpl = TINKER_MOD(potent, use_chgdpl);
int& use_dipole = TINKER_MOD(potent, use_dipole);
int& use_mpole = TINKER_MOD(potent, use_mpole);
int& use_polar = TINKER_MOD(potent, use_polar);
int& use_chgtrn = TINKER_MOD(potent, use_chgtrn);
int& use_chgflx = TINKER_MOD(potent, use_chgflx);
int& use_rxnfld = TINKER_MOD(potent, use_rxnfld);
int& use_solv = TINKER_MOD(potent, use_solv);
int& use_metal = TINKER_MOD(potent, use_metal);
int& use_geom = TINKER_MOD(potent, use_geom);
int& use_extra = TINKER_MOD(potent, use_extra);
int& use_born = TINKER_MOD(potent, use_born);
int& use_orbit = TINKER_MOD(potent, use_orbit);
int& use_mutate = TINKER_MOD(potent, use_mutate);
#endif
} }
