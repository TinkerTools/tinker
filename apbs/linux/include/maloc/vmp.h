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
 * rcsid="$Id: vmp.h,v 1.10 2006/06/03 07:22:30 mholst Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     vmp.h    < vmp.c >
 *
 * Purpose:  Class Vmp: a Virtual MPI communication layer object.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#ifndef _VMP_H_
#define _VMP_H_

#include <maloc/maloc_base.h>

#include <maloc/vsys.h>
#include <maloc/vmpi.h>
#include <maloc/vcom.h>

/*
 * ***************************************************************************
 * Class Vmp: Parameters and datatypes
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * Class Vmp: Definition
 * ***************************************************************************
 */

typedef struct Vmp {
    int  mpi_rank;     /* my process ID                                     */
    int  mpi_size;     /* number of processess in this execution            */
} Vmp;

/*
 * ***************************************************************************
 * Class Vmp: Inlineable methods (vmp.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_MALOC)
#else /* if defined(VINLINE_MALOC) */
#endif /* if !defined(VINLINE_MALOC) */

/*
 * ***************************************************************************
 * Class Vmp: Non-inlineable methods (vmp.c)
 * ***************************************************************************
 */

VEXTERNC int Vmp_init(int *argc, char ***argv);
VEXTERNC int Vmp_finalize(void);

VEXTERNC Vmp* Vmp_ctor(void);
VEXTERNC void Vmp_dtor(Vmp **thee);

VEXTERNC int Vmp_rank(Vmp *thee);
VEXTERNC int Vmp_size(Vmp *thee);
VEXTERNC int Vmp_barr(Vmp *thee);

VEXTERNC int Vmp_send(Vmp *thee, int des, char *buf, int bufsize);
VEXTERNC int Vmp_recv(Vmp *thee, int src, char *buf, int bufsize);

#endif /* _VMP_H_ */

