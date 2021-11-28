/**
 *  @file       psh.h
 *  @brief      Header file for a simple parallel extension of ALOC's VSH.
 *  @version    $Id: psh.h,v 1.21 2008/02/05 00:10:01 fetk Exp $
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

