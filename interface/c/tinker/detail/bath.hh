#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxnose 4
extern int TINKER_MOD(bath, voltrial);
extern double TINKER_MOD(bath, kelvin);
extern double TINKER_MOD(bath, atmsph);
extern double TINKER_MOD(bath, tautemp);
extern double TINKER_MOD(bath, taupres);
extern double TINKER_MOD(bath, compress);
extern double TINKER_MOD(bath, collide);
extern double TINKER_MOD(bath, eta);
extern double TINKER_MOD(bath, volmove);
extern double TINKER_MOD(bath, vbar);
extern double TINKER_MOD(bath, qbar);
extern double TINKER_MOD(bath, gbar);
extern double TINKER_MOD(bath, vnh)[TINKER_MOD__maxnose];
extern double TINKER_MOD(bath, qnh)[TINKER_MOD__maxnose];
extern double TINKER_MOD(bath, gnh)[TINKER_MOD__maxnose];
extern int TINKER_MOD(bath, isothermal);
extern int TINKER_MOD(bath, isobaric);
extern int TINKER_MOD(bath, anisotrop);
extern char TINKER_MOD(bath, volscale)[9];
extern char TINKER_MOD(bath, barostat)[11];
extern char TINKER_MOD(bath, thermostat)[11];
#ifdef __cplusplus
}
#endif
