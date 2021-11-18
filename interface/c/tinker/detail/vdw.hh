#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(vdw, nvdw);
extern int* TINKER_MOD(vdw, ivdw);
extern int* TINKER_MOD(vdw, jvdw);
extern int* TINKER_MOD(vdw, ired);
extern double* TINKER_MOD(vdw, kred);
extern double* TINKER_MOD(vdw, xred);
extern double* TINKER_MOD(vdw, yred);
extern double* TINKER_MOD(vdw, zred);
extern double* TINKER_MOD(vdw, radmin);
extern double* TINKER_MOD(vdw, epsilon);
extern double* TINKER_MOD(vdw, radmin4);
extern double* TINKER_MOD(vdw, epsilon4);
extern double* TINKER_MOD(vdw, radhbnd);
extern double* TINKER_MOD(vdw, epshbnd);
#ifdef __cplusplus
}
#endif
