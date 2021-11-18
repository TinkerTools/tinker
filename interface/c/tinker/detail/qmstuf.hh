#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(qmstuf, ngatom);
extern double TINKER_MOD(qmstuf, egau);
extern double* TINKER_MOD(qmstuf, gx);
extern double* TINKER_MOD(qmstuf, gy);
extern double* TINKER_MOD(qmstuf, gz);
extern double* TINKER_MOD(qmstuf, gfreq);
extern double* TINKER_MOD(qmstuf, gforce);
extern double* TINKER_MOD(qmstuf, gh);
#ifdef __cplusplus
}
#endif
