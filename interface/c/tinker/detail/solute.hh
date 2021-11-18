#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern double TINKER_MOD(solute, doffset);
extern double TINKER_MOD(solute, onipr);
extern double TINKER_MOD(solute, p1);
extern double TINKER_MOD(solute, p2);
extern double TINKER_MOD(solute, p3);
extern double TINKER_MOD(solute, p4);
extern double TINKER_MOD(solute, p5);
extern double* TINKER_MOD(solute, rsolv);
extern double* TINKER_MOD(solute, rdescr);
extern double* TINKER_MOD(solute, asolv);
extern double* TINKER_MOD(solute, rborn);
extern double* TINKER_MOD(solute, drb);
extern double* TINKER_MOD(solute, drbp);
extern double* TINKER_MOD(solute, drobc);
extern double* TINKER_MOD(solute, gpol);
extern double* TINKER_MOD(solute, shct);
extern double* TINKER_MOD(solute, aobc);
extern double* TINKER_MOD(solute, bobc);
extern double* TINKER_MOD(solute, gobc);
extern double* TINKER_MOD(solute, vsolv);
extern double* TINKER_MOD(solute, wace);
extern double* TINKER_MOD(solute, s2ace);
extern double* TINKER_MOD(solute, uace);
#ifdef __cplusplus
}
#endif
