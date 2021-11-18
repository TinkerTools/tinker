#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxpss 500
extern int TINKER_MOD(tree, nlevel);
extern double TINKER_MOD(tree, etree);
extern double TINKER_MOD(tree, ilevel)[TINKER_MOD__maxpss+1];
#ifdef __cplusplus
}
#endif
