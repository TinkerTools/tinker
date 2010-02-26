/*
 * ***************************************************************************
 * MALOC = < Minimal Abstraction Layer for Object-oriented C >
 * Copyright (C) 1994--2000  Michael Holst
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 675 Mass Ave, Cambridge, MA 02139, USA.
 * 
 * rcsid="$Id: vnm.h,v 1.12 2002/10/01 21:29:45 mholst Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     vnm.h    < vnm.c >
 *
 * Purpose:  Header file for an ISO C [V]irtual [N]umerical [M]achine.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#ifndef _VNM_H_
#define _VNM_H_

#include <maloc/maloc_base.h>

/*
 * ***************************************************************************
 * Useful constants and functions (timers, epsilon, token generators, i/o)
 * ***************************************************************************
 */

VEXTERNC int Vnm_sigInt(void);
VEXTERNC void Vnm_sigIntSet(void);
VEXTERNC void Vnm_sigIntClear(void);

VEXTERNC int Vnm_jmpOk(void);
VEXTERNC void Vnm_jmpOkSet(void);
VEXTERNC void Vnm_jmpOkClear(void);

VEXTERNC jmp_buf *Vnm_signalInit(void);
VEXTERNC void Vnm_regHand(void);
VEXTERNC void Vnm_sigHand(int num);

#define VPOW_SAFE(x,y) (Vnm_powsafe(x,y))
VEXTERNC double Vnm_powsafe(double x, double y);
VEXTERNC void Vnm_typeChk(void);
VEXTERNC double Vnm_epsmac(void);
VEXTERNC int Vnm_gentokens(char *buf, char **argv, 
    const int argvmax, const char *white, const char *comment);

#define VTIMERS 100
VEXTERNC void Vnm_tstart(int timer, const char *name);
VEXTERNC void Vnm_tstop(int timer, const char *name);

VEXTERNC char *Vnm_getuser(char *user, int usermax);
VEXTERNC char *Vnm_getos(char *os, int osmax);
VEXTERNC char *Vnm_gethost(char *host, int hostmax);
VEXTERNC char *Vnm_gethome(char *path, int pathmax);
VEXTERNC char *Vnm_getcwd(char *path, int pathmax);
VEXTERNC int Vnm_chdir(const char *path);
VEXTERNC int Vnm_mkdir(const char *path);

VEXTERNC int Vnm_system(const char *cmd);
VEXTERNC int Vnm_systemBack(const char *cmd);
VEXTERNC int Vnm_systemKill(const char *cmd);
VEXTERNC int Vnm_exec(int argc, char **argv);
VEXTERNC void Vnm_sleep(int nusecs);

VEXTERNC int Vnm_ioTag(void);
VEXTERNC int Vnm_nTags(void);
VEXTERNC void Vnm_setIoTag(int myTag, int numTags);
VEXTERNC FILE *Vnm_open(const int unit);
VEXTERNC int Vnm_close(const int unit);
VEXTERNC void Vnm_flush(const int unit);
VEXTERNC void Vnm_redirect(const int flag);
VEXTERNC void Vnm_print(const int unit, const char *format, ...);
VEXTERNC void Vnm_tprint(const int unit, const char *format, ...);

VEXTERNC void Vnm_qsort(int *u, int size);
VEXTERNC void Vnm_qsortOrd(int *u, int *ord, int size);
VEXTERNC void Vnm_dqsort(double *u, int size);
VEXTERNC void Vnm_dqsortOrd(double *u, int *ord, int size);

#endif /* _VNM_H_ */

