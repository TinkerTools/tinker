#pragma once

#include "macro.hh"

namespace tinker { namespace polgrp {
const int maxp11 = 200;
const int maxp12 = 200;
const int maxp13 = 200;
const int maxp14 = 200;
extern int*& np11;
extern int*& np12;
extern int*& np13;
extern int*& np14;
extern int*& ip11;
extern int*& ip12;
extern int*& ip13;
extern int*& ip14;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int* TINKER_MOD(polgrp, np11);
extern "C" int* TINKER_MOD(polgrp, np12);
extern "C" int* TINKER_MOD(polgrp, np13);
extern "C" int* TINKER_MOD(polgrp, np14);
extern "C" int* TINKER_MOD(polgrp, ip11);
extern "C" int* TINKER_MOD(polgrp, ip12);
extern "C" int* TINKER_MOD(polgrp, ip13);
extern "C" int* TINKER_MOD(polgrp, ip14);

int*& np11 = TINKER_MOD(polgrp, np11);
int*& np12 = TINKER_MOD(polgrp, np12);
int*& np13 = TINKER_MOD(polgrp, np13);
int*& np14 = TINKER_MOD(polgrp, np14);
int*& ip11 = TINKER_MOD(polgrp, ip11);
int*& ip12 = TINKER_MOD(polgrp, ip12);
int*& ip13 = TINKER_MOD(polgrp, ip13);
int*& ip14 = TINKER_MOD(polgrp, ip14);
#endif
} }
