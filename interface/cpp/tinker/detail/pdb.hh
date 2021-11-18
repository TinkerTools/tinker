#pragma once

#include "macro.hh"

namespace tinker { namespace pdb {
extern int& npdb;
extern int& nres;
extern int*& resnum;
extern int*& resatm;
extern int*& npdb12;
extern int*& ipdb12;
extern int*& pdblist;
extern double*& xpdb;
extern double*& ypdb;
extern double*& zpdb;
extern char (&altsym)[1];
extern char (*&pdbres)[3];
extern char (*&pdbsym)[3];
extern char (*&pdbatm)[4];
extern char (*&pdbtyp)[6];
extern char (&chnsym)[20];
extern char (&instyp)[20];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(pdb, npdb);
extern "C" int TINKER_MOD(pdb, nres);
extern "C" int* TINKER_MOD(pdb, resnum);
extern "C" int* TINKER_MOD(pdb, resatm);
extern "C" int* TINKER_MOD(pdb, npdb12);
extern "C" int* TINKER_MOD(pdb, ipdb12);
extern "C" int* TINKER_MOD(pdb, pdblist);
extern "C" double* TINKER_MOD(pdb, xpdb);
extern "C" double* TINKER_MOD(pdb, ypdb);
extern "C" double* TINKER_MOD(pdb, zpdb);
extern "C" char TINKER_MOD(pdb, altsym)[1];
extern "C" char (*TINKER_MOD(pdb, pdbres))[3];
extern "C" char (*TINKER_MOD(pdb, pdbsym))[3];
extern "C" char (*TINKER_MOD(pdb, pdbatm))[4];
extern "C" char (*TINKER_MOD(pdb, pdbtyp))[6];
extern "C" char TINKER_MOD(pdb, chnsym)[20];
extern "C" char TINKER_MOD(pdb, instyp)[20];

int& npdb = TINKER_MOD(pdb, npdb);
int& nres = TINKER_MOD(pdb, nres);
int*& resnum = TINKER_MOD(pdb, resnum);
int*& resatm = TINKER_MOD(pdb, resatm);
int*& npdb12 = TINKER_MOD(pdb, npdb12);
int*& ipdb12 = TINKER_MOD(pdb, ipdb12);
int*& pdblist = TINKER_MOD(pdb, pdblist);
double*& xpdb = TINKER_MOD(pdb, xpdb);
double*& ypdb = TINKER_MOD(pdb, ypdb);
double*& zpdb = TINKER_MOD(pdb, zpdb);
char (&altsym)[1] = TINKER_MOD(pdb, altsym);
char (*&pdbres)[3] = TINKER_MOD(pdb, pdbres);
char (*&pdbsym)[3] = TINKER_MOD(pdb, pdbsym);
char (*&pdbatm)[4] = TINKER_MOD(pdb, pdbatm);
char (*&pdbtyp)[6] = TINKER_MOD(pdb, pdbtyp);
char (&chnsym)[20] = TINKER_MOD(pdb, chnsym);
char (&instyp)[20] = TINKER_MOD(pdb, instyp);
#endif
} }
