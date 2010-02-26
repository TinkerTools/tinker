/*
 * ***************************************************************************
 * MALOC = < Minimal Abstraction Layer for Object-oriented C >
 * Copyright (C) 1994--2006  Michael Holst
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
 * rcsid="$Id: psh.h,v 1.19 2006/06/03 07:22:30 mholst Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     psh.h    < psh.c, ... >
 *
 * Purpose:  Header file for a simple parallel extension of ALOC's VSH.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#ifndef _PSH_H_
#define _PSH_H_

#include <maloc/maloc_base.h>

#include <maloc/vsys.h>
#include <maloc/vsh.h>
#include <maloc/vmp.h>

VEXTERNC int Vsh_pshell(Vsh *thee, char *pPR, void *pthee,
    int (*builtin)(void *thee, int argc, char **argv));

VEXTERNC Vio *Vsh_pioSetup(Vsh *thee, char *key);
VEXTERNC void Vsh_pioCleanup(Vsh *thee, Vio **sock);

#endif /* _PSH_H_ */

