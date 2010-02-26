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
 * rcsid="$Id: vmpi.h,v 1.18 2006/06/03 07:22:30 mholst Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     vmpi.h    < vmpi.c >
 *
 * Purpose:  Class Vmpi: a Virtual MPI communication layer object.
 *
 * Notes:    Class Vmpi is a thin object-oriented Clean C layer on top of the
 *           MPI communication library.  Vmpi provides access to the minimal
 *           set of ten MPI primitives required to implement the Bank-Holst
 *           parallel adaptive algorithm, using either the Bank-Holst Oracle
 *           library, or directly.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#ifndef _VMPI_H_
#define _VMPI_H_

#include <maloc/maloc_base.h>

#include <maloc/vsys.h>

/*
 * ***************************************************************************
 * Class Vmpi: Parameters and datatypes
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * Class Vmpi: Definition
 * ***************************************************************************
 */

typedef struct Vmpi {
    int  mpi_rank;     /* my process ID                                     */
    int  mpi_size;     /* number of processess in this execution            */
} Vmpi;

/*
 * ***************************************************************************
 * Class Vmpi: Inlineable methods (vmpi.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_MALOC)
#else /* if defined(VINLINE_MALOC) */
#endif /* if !defined(VINLINE_MALOC) */

/*
 * ***************************************************************************
 * Class Vmpi: Non-inlineable methods (vmpi.c)
 * ***************************************************************************
 */

VEXTERNC int Vmpi_init(int *argc, char ***argv);
VEXTERNC int Vmpi_finalize(void);

VEXTERNC Vmpi* Vmpi_ctor(void);
VEXTERNC void Vmpi_dtor(Vmpi **thee);

VEXTERNC int Vmpi_rank(Vmpi *thee);
VEXTERNC int Vmpi_size(Vmpi *thee);
VEXTERNC int Vmpi_barr(Vmpi *thee);

VEXTERNC int Vmpi_send(Vmpi *thee, int des, char *buf, int bufsize);
VEXTERNC int Vmpi_recv(Vmpi *thee, int src, char *buf, int bufsize);

VEXTERNC int Vmpi_bcast(Vmpi *thee, char *buf, int bufsize);
VEXTERNC int Vmpi_reduce(Vmpi *thee, char *sbuf, char *rbuf, int bufsize);
VEXTERNC int Vmpi_isend(Vmpi *thee, int des, char *buf, int bufsize);

#endif /* _VMPI_H_ */

