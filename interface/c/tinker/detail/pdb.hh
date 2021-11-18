#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(pdb, npdb);
extern int TINKER_MOD(pdb, nres);
extern int* TINKER_MOD(pdb, resnum);
extern int* TINKER_MOD(pdb, resatm);
extern int* TINKER_MOD(pdb, npdb12);
extern int* TINKER_MOD(pdb, ipdb12);
extern int* TINKER_MOD(pdb, pdblist);
extern double* TINKER_MOD(pdb, xpdb);
extern double* TINKER_MOD(pdb, ypdb);
extern double* TINKER_MOD(pdb, zpdb);
extern char TINKER_MOD(pdb, altsym)[1];
extern char (*TINKER_MOD(pdb, pdbres))[3];
extern char (*TINKER_MOD(pdb, pdbsym))[3];
extern char (*TINKER_MOD(pdb, pdbatm))[4];
extern char (*TINKER_MOD(pdb, pdbtyp))[6];
extern char TINKER_MOD(pdb, chnsym)[20];
extern char TINKER_MOD(pdb, instyp)[20];
#ifdef __cplusplus
}
#endif
