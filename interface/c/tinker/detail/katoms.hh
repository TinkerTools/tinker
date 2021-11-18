#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int* TINKER_MOD(katoms, atmcls);
extern int* TINKER_MOD(katoms, atmnum);
extern int* TINKER_MOD(katoms, ligand);
extern double* TINKER_MOD(katoms, weight);
extern char (*TINKER_MOD(katoms, symbol))[3];
extern char (*TINKER_MOD(katoms, describe))[24];
#ifdef __cplusplus
}
#endif
