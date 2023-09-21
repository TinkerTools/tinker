#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__epso 0.1100e0
#define TINKER_MOD__epsh 0.0135e0
#define TINKER_MOD__rmino 1.7025e0
#define TINKER_MOD__rminh 1.3275e0
#define TINKER_MOD__awater 0.033428e0
#define TINKER_MOD__slevy 1.0e0
#define TINKER_MOD__shctd 0.75e0
#define TINKER_MOD__cavoff 0.0e0
#define TINKER_MOD__dspoff 1.056e0
extern double TINKER_MOD(nonpol, solvprs);
extern double TINKER_MOD(nonpol, surften);
extern double TINKER_MOD(nonpol, spcut);
extern double TINKER_MOD(nonpol, spoff);
extern double TINKER_MOD(nonpol, stcut);
extern double TINKER_MOD(nonpol, stoff);
extern double* TINKER_MOD(nonpol, radcav);
extern double* TINKER_MOD(nonpol, raddsp);
extern double* TINKER_MOD(nonpol, epsdsp);
extern double* TINKER_MOD(nonpol, cdsp);
#ifdef __cplusplus
}
#endif
