#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern double TINKER_MOD(moment, netchg);
extern double TINKER_MOD(moment, netdpl);
extern double TINKER_MOD(moment, netqpl)[3];
extern double TINKER_MOD(moment, xdpl);
extern double TINKER_MOD(moment, ydpl);
extern double TINKER_MOD(moment, zdpl);
extern double TINKER_MOD(moment, xxqpl);
extern double TINKER_MOD(moment, xyqpl);
extern double TINKER_MOD(moment, xzqpl);
extern double TINKER_MOD(moment, yxqpl);
extern double TINKER_MOD(moment, yyqpl);
extern double TINKER_MOD(moment, yzqpl);
extern double TINKER_MOD(moment, zxqpl);
extern double TINKER_MOD(moment, zyqpl);
extern double TINKER_MOD(moment, zzqpl);
extern int* TINKER_MOD(moment, momuse);
#ifdef __cplusplus
}
#endif
