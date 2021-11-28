/**
 *  @file       vsh.h
 *  @brief      Header file for vsh, a bourne-compatible shell.
 *  @version    $Id: vsh.h,v 1.20 2008/03/12 05:13:59 fetk Exp $
 *  @author     Michael Holst
 *
 *  @attention
 *  @verbatim
 *

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


#ifndef _VSH_H_
#define _VSH_H_

#include <maloc/maloc_base.h>

#include <maloc/vsys.h>

/*
 * ***************************************************************************
 * Class Vsh: Parameters and datatypes
 * ***************************************************************************
 */

/** @brief Class Vsh: Definition */

typedef struct Vsh {

    Vmem   *vmem;        /**< the memory manager                            */
    int    iMadeVmem;    /**< did i make vmem or was it inherited           */

    char processArgs;    /**< whether the shell should process (argc,argv)  */

    int envValuLen;      /**< number of environment variables               */
    int envInfoLen;      /**< number of environment variable help strings   */
    char **envValu;      /**< the environment variables                     */
    char **envInfo;      /**< the environment variable help strings         */

    FILE *inUnit;        /**< input unit                                    */
    FILE *scUnit;        /**< script input unit                             */
    FILE *clUnit;        /**< input unit                                    */
    FILE *cinUnit;       /**< input unit                                    */
    char cinName[VMAX_ARGLEN];  /**< input unit                             */

    char PR[VMAX_ARGLEN];       /**< minimal prompt (just the binary name)  */
    char PR_PATH[VMAX_ARGLEN];  /**< full prompt (user,hostname,path,etc)   */
    char PR_EXIT[VMAX_ARGLEN];  /**< the exit print string                  */

    int cmdKey;          /**< external supershell command key               */
    void *Ext_thee;      /**< external supershell object                    */

    char *buf;           /**< internal buffer                               */
    int bufsize;         /**< internal buffer size                          */

    /** @brief external supershell builtin function */
    int (*Ext_builtin)(void *thee, int argc, char **argv);

} Vsh;                                                                        

/*
 * ***************************************************************************
 * Class Vsh: Inlineable methods (vsh.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_MALOC)
#else /* if defined(VINLINE_MALOC) */
#endif /* if !defined(VINLINE_MALOC) */


/** @brief Class Vsh: Non-inlineable method (vsh.c) */
VEXTERNC Vsh* Vsh_ctor(Vmem *vmem, int argc, char **argv);
/** @brief Class Vsh: Non-inlineable method (vsh.c) */
VEXTERNC void Vsh_dtor(Vsh **thee);

/** @brief Class Vsh: Non-inlineable method (vsh.c) */
VEXTERNC int Vsh_shell(Vsh *thee, char *pPR, void *pthee,
    int (*builtin)(void *thee, int argc, char **argv));

/** @brief Class Vsh: Non-inlineable method (vsh.c) */
VEXTERNC int Vsh_putenv(Vsh *thee, const char *envi, const char *valu);
/** @brief Class Vsh: Non-inlineable method (vsh.c) */
VEXTERNC int Vsh_putenvInfo(Vsh *thee, const char *envi, const char *valu);
/** @brief Class Vsh: Non-inlineable method (vsh.c) */
VEXTERNC int Vsh_putenvInt(Vsh *thee, const char *envi, const int valu);
/** @brief Class Vsh: Non-inlineable method (vsh.c) */
VEXTERNC int Vsh_putenvReal(Vsh *thee, const char *envi, const double valu);

/** @brief Class Vsh: Non-inlineable method (vsh.c) */
VEXTERNC char *Vsh_getenv(Vsh *thee, const char *envi);
/** @brief Class Vsh: Non-inlineable method (vsh.c) */
VEXTERNC char *Vsh_getenvInfo(Vsh *thee, const char *envi);
/** @brief Class Vsh: Non-inlineable method (vsh.c) */
VEXTERNC int Vsh_getenvInt(Vsh *thee, const char *envi);
/** @brief Class Vsh: Non-inlineable method (vsh.c) */
VEXTERNC double Vsh_getenvReal(Vsh *thee, const char *envi);

/** @brief Class Vsh: Non-inlineable method (vsh.c) */
VEXTERNC void Vsh_remove(Vsh *thee, const char *envi);
/** @brief Class Vsh: Non-inlineable method (vsh.c) */
VEXTERNC void Vsh_wipe(Vsh *thee);

/** @brief Class Vsh: Non-inlineable method (vsh.c) */
VEXTERNC void Vsh_memChk(Vsh *thee);

/** @brief Class Vsh: Non-inlineable method (vsh.c) */
VEXTERNC Vio *Vsh_ioSetup(Vsh *thee, char *key);
/** @brief Class Vsh: Non-inlineable method (vsh.c) */
VEXTERNC void Vsh_ioCleanup(Vsh *thee, Vio **sock);

#endif /* _VSH_H_ */

