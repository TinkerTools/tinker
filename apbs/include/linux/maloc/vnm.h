/**
 *  @file       vnm.h
 *  @brief      Header file for an ISO C [V]irtual [N]umerical [M]achine.
 *  @version    $Id: vnm.h,v 1.15 2008/02/05 00:11:34 fetk Exp $
 *  @author     Michael Holst
 *  
 *  @attention
 *  @verbatim
 *
 * MALOC = < Minimal Abstraction Layer for Object-oriented C >
 * Copyright (C) 1994--2008 Michael Holst
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * 
 *  @endverbatim
 */


#ifndef _VNM_H_
#define _VNM_H_

#include <maloc/maloc_base.h>


/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC int Vnm_sigInt(void);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC void Vnm_sigIntSet(void);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC void Vnm_sigIntClear(void);

/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC int Vnm_jmpOk(void);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC void Vnm_jmpOkSet(void);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC void Vnm_jmpOkClear(void);

/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC jmp_buf *Vnm_signalInit(void);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC void Vnm_regHand(void);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC void Vnm_sigHand(int num);

/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
#define VPOW_SAFE(x,y) (Vnm_powsafe(x,y))
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC double Vnm_powsafe(double x, double y);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC void Vnm_typeChk(void);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC double Vnm_epsmac(void);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC int Vnm_gentokens(char *buf, char **argv, 
    const int argvmax, const char *white, const char *comment);

/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
#define VTIMERS 100
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC void Vnm_tstart(int timer, const char *name);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC void Vnm_tstop(int timer, const char *name);

/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC char *Vnm_getuser(char *user, int usermax);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC char *Vnm_getos(char *os, int osmax);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC char *Vnm_gethost(char *host, int hostmax);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC char *Vnm_gethome(char *path, int pathmax);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC char *Vnm_getcwd(char *path, int pathmax);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC int Vnm_chdir(const char *path);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC int Vnm_mkdir(const char *path);

/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC int Vnm_system(const char *cmd);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC int Vnm_systemBack(const char *cmd);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC int Vnm_systemKill(const char *cmd);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC int Vnm_exec(int argc, char **argv);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC void Vnm_sleep(int nusecs);

/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC int Vnm_ioTag(void);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC int Vnm_nTags(void);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC void Vnm_setIoTag(int myTag, int numTags);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC FILE *Vnm_open(const int unit);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC int Vnm_close(const int unit);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC void Vnm_flush(const int unit);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC void Vnm_redirect(const int flag);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC void Vnm_print(const int unit, const char *format, ...);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC void Vnm_tprint(const int unit, const char *format, ...);

/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC void Vnm_qsort(int *u, int size);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC void Vnm_qsortOrd(int *u, int *ord, int size);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC void Vnm_dqsort(double *u, int size);
/** @brief Useful constants and functions (timers, epsilon, token generators, i/o) */
VEXTERNC void Vnm_dqsortOrd(double *u, int *ord, int size);

#endif /* _VNM_H_ */

