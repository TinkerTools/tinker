#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(minima, maxiter);
extern int TINKER_MOD(minima, nextiter);
extern double TINKER_MOD(minima, fctmin);
extern double TINKER_MOD(minima, hguess);
#ifdef __cplusplus
}
#endif
