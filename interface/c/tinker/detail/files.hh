#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(files, nprior);
extern int TINKER_MOD(files, ldir);
extern int TINKER_MOD(files, leng);
extern char TINKER_MOD(files, filename)[240];
extern char TINKER_MOD(files, outfile)[240];
#ifdef __cplusplus
}
#endif
