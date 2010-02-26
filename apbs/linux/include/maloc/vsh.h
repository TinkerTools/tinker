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
 * rcsid="$Id: vsh.h,v 1.15 2002/10/01 21:29:45 mholst Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     vsh.h    < vsh.c, ... >
 *
 * Purpose:  Header file for vsh, a bourne-compatible shell.
 *
 * Author:   Michael Holst
 * ***************************************************************************
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

/*
 * ***************************************************************************
 * Class Vsh: Definition
 * ***************************************************************************
 */

typedef struct Vsh {

    Vmem   *vmem;        /* the memory manager                              */
    int    iMadeVmem;    /* did i make vmem or was it inherited             */

    char processArgs;    /* whether the shell should process (argc,argv)    */

    int envValuLen;      /* number of environment variables                 */
    int envInfoLen;      /* number of environment variable help strings     */
    char **envValu;      /* the environment variables                       */
    char **envInfo;      /* the environment variable help strings           */

    FILE *inUnit;        /* input unit                                      */
    FILE *scUnit;        /* script input unit                               */
    FILE *clUnit;        /* input unit                                      */
    FILE *cinUnit;       /* input unit                                      */
    char cinName[80];    /* input unit                                      */

    char PR[80];         /* minimal prompt (just the binary name)           */
    char PR_PATH[80];    /* full prompt (with user,hostname,path,etc)       */
    char PR_EXIT[80];    /* the exit print string                           */

    int cmdKey;          /* external supershell command key                 */
    void *Ext_thee;      /* external supershell object                      */

    char *buf;           /* internal buffer                                 */
    int bufsize;         /* internal buffer size                            */

    /* external supershell builtin function */
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

/*
 * ***************************************************************************
 * Class Vsh: Non-inlineable methods (vsh.c)
 * ***************************************************************************
 */

VEXTERNC Vsh* Vsh_ctor(Vmem *vmem, int argc, char **argv);
VEXTERNC void Vsh_dtor(Vsh **thee);

VEXTERNC int Vsh_shell(Vsh *thee, char *pPR, void *pthee,
    int (*builtin)(void *thee, int argc, char **argv));

VEXTERNC int Vsh_putenv(Vsh *thee, const char *envi, const char *valu);
VEXTERNC int Vsh_putenvInfo(Vsh *thee, const char *envi, const char *valu);
VEXTERNC int Vsh_putenvInt(Vsh *thee, const char *envi, const int valu);
VEXTERNC int Vsh_putenvReal(Vsh *thee, const char *envi, const double valu);

VEXTERNC char *Vsh_getenv(Vsh *thee, const char *envi);
VEXTERNC char *Vsh_getenvInfo(Vsh *thee, const char *envi);
VEXTERNC int Vsh_getenvInt(Vsh *thee, const char *envi);
VEXTERNC double Vsh_getenvReal(Vsh *thee, const char *envi);

VEXTERNC void Vsh_remove(Vsh *thee, const char *envi);
VEXTERNC void Vsh_wipe(Vsh *thee);

VEXTERNC void Vsh_memChk(Vsh *thee);

VEXTERNC Vio *Vsh_ioSetup(Vsh *thee, char *key);
VEXTERNC void Vsh_ioCleanup(Vsh *thee, Vio **sock);

#endif /* _VSH_H_ */

