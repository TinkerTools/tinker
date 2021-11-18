#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(dipole, ndipole);
extern int* TINKER_MOD(dipole, idpl);
extern double* TINKER_MOD(dipole, bdpl);
extern double* TINKER_MOD(dipole, sdpl);
#ifdef __cplusplus
}
#endif
