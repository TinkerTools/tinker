#pragma once

#include "macro.hh"

namespace tinker { namespace light {
extern int& nlight;
extern int*& kbx;
extern int*& kby;
extern int*& kbz;
extern int*& kex;
extern int*& key;
extern int*& kez;
extern int*& locx;
extern int*& locy;
extern int*& locz;
extern int*& rgx;
extern int*& rgy;
extern int*& rgz;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(light, nlight);
extern "C" int* TINKER_MOD(light, kbx);
extern "C" int* TINKER_MOD(light, kby);
extern "C" int* TINKER_MOD(light, kbz);
extern "C" int* TINKER_MOD(light, kex);
extern "C" int* TINKER_MOD(light, key);
extern "C" int* TINKER_MOD(light, kez);
extern "C" int* TINKER_MOD(light, locx);
extern "C" int* TINKER_MOD(light, locy);
extern "C" int* TINKER_MOD(light, locz);
extern "C" int* TINKER_MOD(light, rgx);
extern "C" int* TINKER_MOD(light, rgy);
extern "C" int* TINKER_MOD(light, rgz);

int& nlight = TINKER_MOD(light, nlight);
int*& kbx = TINKER_MOD(light, kbx);
int*& kby = TINKER_MOD(light, kby);
int*& kbz = TINKER_MOD(light, kbz);
int*& kex = TINKER_MOD(light, kex);
int*& key = TINKER_MOD(light, key);
int*& kez = TINKER_MOD(light, kez);
int*& locx = TINKER_MOD(light, locx);
int*& locy = TINKER_MOD(light, locy);
int*& locz = TINKER_MOD(light, locz);
int*& rgx = TINKER_MOD(light, rgx);
int*& rgy = TINKER_MOD(light, rgy);
int*& rgz = TINKER_MOD(light, rgz);
#endif
} }
