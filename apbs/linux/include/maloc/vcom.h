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
 * rcsid="$Id: vcom.h,v 1.26 2006/06/03 07:22:30 mholst Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     vcom.h    < vcom.c >
 *
 * Purpose:  Class Vcom: virtual (currently just MPI) communications layer
 *
 * Author:   Nathan Baker and Michael Holst
 * ***************************************************************************
 */

#ifndef _VCOM_H_
#define _VCOM_H_

#include <maloc/maloc_base.h>

#include <maloc/vsys.h>

/* A base value for MPI tags */
#define VCOM_MPI_TAG 111

/*
 * ***************************************************************************
 * Class Vcom: Parameters and datatypes
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * Class Vcom: Definition
 * ***************************************************************************
 */

typedef struct Vcom {

    int  mpi_rank;     /* Local PE rank from MPI                            */
    int  mpi_size;     /* Total number of PEs in this communicator from MPI */ 

    int  type;         /* Communications type                               */
                       /*   0 = not initialized                             */
                       /*   1 = Message Passing Interface 1.1               */

    int  error;        /* note if any error has occurred on this vcom device*/

    void *core;        /* Private MPI core                                  */

} Vcom;

/*
 * ***************************************************************************
 * Class Vcom: Inlineable methods (vcom.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_MALOC)
#else /* if defined(VINLINE_MALOC) */
#endif /* if !defined(VINLINE_MALOC) */

/*
 * ***************************************************************************
 * Class Vcom: Non-Inlineable methods (vcom.c)
 * ***************************************************************************
 */

VEXTERNC int Vcom_init(int *argc, char ***argv);
VEXTERNC int Vcom_finalize(void);

VEXTERNC Vcom* Vcom_ctor(int commtype);
VEXTERNC int Vcom_ctor2(Vcom* thee, int commtype);
VEXTERNC void Vcom_dtor(Vcom **thee);
VEXTERNC void Vcom_dtor2(Vcom *thee);

VEXTERNC int Vcom_send(Vcom *thee, int des, void *buf, int len, int type, 
  int block);
VEXTERNC int Vcom_recv(Vcom *thee, int src, void *buf, int len, int type, 
  int block);
VEXTERNC int Vcom_getCount(Vcom *thee, int src, int *length, int type);
VEXTERNC int Vcom_reduce(Vcom *thee, void *sendbuf, void *recvbuf, int length, 
  int type, int op);
VEXTERNC int Vcom_size(Vcom *thee);
VEXTERNC int Vcom_resize(Vcom *thee, int newsize);
VEXTERNC int Vcom_rank(Vcom *thee);
VEXTERNC int Vcom_barr(Vcom *thee);

#endif /* _VCOM_H_ */

