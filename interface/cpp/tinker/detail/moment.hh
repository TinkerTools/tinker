#pragma once

#include "macro.hh"

namespace tinker { namespace moment {
extern double& netchg;
extern double& netdpl;
extern double (&netqpl)[3];
extern double& xdpl;
extern double& ydpl;
extern double& zdpl;
extern double& xxqpl;
extern double& xyqpl;
extern double& xzqpl;
extern double& yxqpl;
extern double& yyqpl;
extern double& yzqpl;
extern double& zxqpl;
extern double& zyqpl;
extern double& zzqpl;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(moment, netchg);
extern "C" double TINKER_MOD(moment, netdpl);
extern "C" double TINKER_MOD(moment, netqpl)[3];
extern "C" double TINKER_MOD(moment, xdpl);
extern "C" double TINKER_MOD(moment, ydpl);
extern "C" double TINKER_MOD(moment, zdpl);
extern "C" double TINKER_MOD(moment, xxqpl);
extern "C" double TINKER_MOD(moment, xyqpl);
extern "C" double TINKER_MOD(moment, xzqpl);
extern "C" double TINKER_MOD(moment, yxqpl);
extern "C" double TINKER_MOD(moment, yyqpl);
extern "C" double TINKER_MOD(moment, yzqpl);
extern "C" double TINKER_MOD(moment, zxqpl);
extern "C" double TINKER_MOD(moment, zyqpl);
extern "C" double TINKER_MOD(moment, zzqpl);

double& netchg = TINKER_MOD(moment, netchg);
double& netdpl = TINKER_MOD(moment, netdpl);
double (&netqpl)[3] = TINKER_MOD(moment, netqpl);
double& xdpl = TINKER_MOD(moment, xdpl);
double& ydpl = TINKER_MOD(moment, ydpl);
double& zdpl = TINKER_MOD(moment, zdpl);
double& xxqpl = TINKER_MOD(moment, xxqpl);
double& xyqpl = TINKER_MOD(moment, xyqpl);
double& xzqpl = TINKER_MOD(moment, xzqpl);
double& yxqpl = TINKER_MOD(moment, yxqpl);
double& yyqpl = TINKER_MOD(moment, yyqpl);
double& yzqpl = TINKER_MOD(moment, yzqpl);
double& zxqpl = TINKER_MOD(moment, zxqpl);
double& zyqpl = TINKER_MOD(moment, zyqpl);
double& zzqpl = TINKER_MOD(moment, zzqpl);
#endif
} }
