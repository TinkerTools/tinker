#pragma once

#include "macro.hh"

namespace tinker { namespace output {
extern int& nonly;
extern int*& ionly;
extern int*& ionlyinv;
extern double*& print3n;
extern int& archive;
extern int& binary;
extern int& noversion;
extern int& overwrite;
extern int& coordsave;
extern int& dynsave;
extern int& cyclesave;
extern int& onlysave;
extern int& arcsave;
extern int& dcdsave;
extern int& velsave;
extern int& frcsave;
extern int& uindsave;
extern int& ustcsave;
extern int& uchgsave;
extern int& usyssave;
extern int& vsyssave;
extern int& udirsave;
extern int& defsave;
extern int& tefsave;
extern char (&coordtype)[9];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(output, nonly);
extern "C" int* TINKER_MOD(output, ionly);
extern "C" int* TINKER_MOD(output, ionlyinv);
extern "C" double* TINKER_MOD(output, print3n);
extern "C" int TINKER_MOD(output, archive);
extern "C" int TINKER_MOD(output, binary);
extern "C" int TINKER_MOD(output, noversion);
extern "C" int TINKER_MOD(output, overwrite);
extern "C" int TINKER_MOD(output, coordsave);
extern "C" int TINKER_MOD(output, dynsave);
extern "C" int TINKER_MOD(output, cyclesave);
extern "C" int TINKER_MOD(output, onlysave);
extern "C" int TINKER_MOD(output, arcsave);
extern "C" int TINKER_MOD(output, dcdsave);
extern "C" int TINKER_MOD(output, velsave);
extern "C" int TINKER_MOD(output, frcsave);
extern "C" int TINKER_MOD(output, uindsave);
extern "C" int TINKER_MOD(output, ustcsave);
extern "C" int TINKER_MOD(output, uchgsave);
extern "C" int TINKER_MOD(output, usyssave);
extern "C" int TINKER_MOD(output, vsyssave);
extern "C" int TINKER_MOD(output, udirsave);
extern "C" int TINKER_MOD(output, defsave);
extern "C" int TINKER_MOD(output, tefsave);
extern "C" char TINKER_MOD(output, coordtype)[9];

int& nonly = TINKER_MOD(output, nonly);
int*& ionly = TINKER_MOD(output, ionly);
int*& ionlyinv = TINKER_MOD(output, ionlyinv);
double*& print3n = TINKER_MOD(output, print3n);
int& archive = TINKER_MOD(output, archive);
int& binary = TINKER_MOD(output, binary);
int& noversion = TINKER_MOD(output, noversion);
int& overwrite = TINKER_MOD(output, overwrite);
int& coordsave = TINKER_MOD(output, coordsave);
int& dynsave = TINKER_MOD(output, dynsave);
int& cyclesave = TINKER_MOD(output, cyclesave);
int& onlysave = TINKER_MOD(output, onlysave);
int& arcsave = TINKER_MOD(output, arcsave);
int& dcdsave = TINKER_MOD(output, dcdsave);
int& velsave = TINKER_MOD(output, velsave);
int& frcsave = TINKER_MOD(output, frcsave);
int& uindsave = TINKER_MOD(output, uindsave);
int& ustcsave = TINKER_MOD(output, ustcsave);
int& uchgsave = TINKER_MOD(output, uchgsave);
int& usyssave = TINKER_MOD(output, usyssave);
int& vsyssave = TINKER_MOD(output, vsyssave);
int& udirsave = TINKER_MOD(output, udirsave);
int& defsave = TINKER_MOD(output, defsave);
int& tefsave = TINKER_MOD(output, tefsave);
char (&coordtype)[9] = TINKER_MOD(output, coordtype);
#endif
} }
