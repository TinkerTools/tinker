#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(socket, skttyp);
extern int TINKER_MOD(socket, cstep);
extern double TINKER_MOD(socket, cdt);
extern double TINKER_MOD(socket, cenergy);
extern int TINKER_MOD(socket, sktstart);
extern int TINKER_MOD(socket, sktstop);
extern int TINKER_MOD(socket, use_socket);
#ifdef __cplusplus
}
#endif
