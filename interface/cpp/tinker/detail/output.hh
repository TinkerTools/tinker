#pragma once

#include "macro.hh"

namespace tinker { namespace output {
extern int& archive;
extern int& binary;
extern int& noversion;
extern int& overwrite;
extern int& cyclesave;
extern int& arcsave;
extern int& dcdsave;
extern int& velsave;
extern int& frcsave;
extern int& uindsave;
extern char (&coordtype)[9];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(output, archive);
extern "C" int TINKER_MOD(output, binary);
extern "C" int TINKER_MOD(output, noversion);
extern "C" int TINKER_MOD(output, overwrite);
extern "C" int TINKER_MOD(output, cyclesave);
extern "C" int TINKER_MOD(output, arcsave);
extern "C" int TINKER_MOD(output, dcdsave);
extern "C" int TINKER_MOD(output, velsave);
extern "C" int TINKER_MOD(output, frcsave);
extern "C" int TINKER_MOD(output, uindsave);
extern "C" char TINKER_MOD(output, coordtype)[9];

int& archive = TINKER_MOD(output, archive);
int& binary = TINKER_MOD(output, binary);
int& noversion = TINKER_MOD(output, noversion);
int& overwrite = TINKER_MOD(output, overwrite);
int& cyclesave = TINKER_MOD(output, cyclesave);
int& arcsave = TINKER_MOD(output, arcsave);
int& dcdsave = TINKER_MOD(output, dcdsave);
int& velsave = TINKER_MOD(output, velsave);
int& frcsave = TINKER_MOD(output, frcsave);
int& uindsave = TINKER_MOD(output, uindsave);
char (&coordtype)[9] = TINKER_MOD(output, coordtype);
#endif
} }
