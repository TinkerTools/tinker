#pragma once

#include "macro.hh"

namespace tinker { namespace socket {
extern int& skttyp;
extern int& cstep;
extern double& cdt;
extern double& cenergy;
extern int& sktstart;
extern int& sktstop;
extern int& use_socket;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(socket, skttyp);
extern "C" int TINKER_MOD(socket, cstep);
extern "C" double TINKER_MOD(socket, cdt);
extern "C" double TINKER_MOD(socket, cenergy);
extern "C" int TINKER_MOD(socket, sktstart);
extern "C" int TINKER_MOD(socket, sktstop);
extern "C" int TINKER_MOD(socket, use_socket);

int& skttyp = TINKER_MOD(socket, skttyp);
int& cstep = TINKER_MOD(socket, cstep);
double& cdt = TINKER_MOD(socket, cdt);
double& cenergy = TINKER_MOD(socket, cenergy);
int& sktstart = TINKER_MOD(socket, sktstart);
int& sktstop = TINKER_MOD(socket, sktstop);
int& use_socket = TINKER_MOD(socket, use_socket);
#endif
} }
