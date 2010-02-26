/**
 *  @file       vmpi.h
 *  @brief      Class Vmpi: a Virtual MPI communication layer object.
 *  @version    $Id: vmpi.h,v 1.20 2008/02/05 00:10:01 fetk Exp $
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


#ifndef _VMPI_H_
#define _VMPI_H_

#include <maloc/maloc_base.h>

#include <maloc/vsys.h>

/*
 * ***************************************************************************
 * Class Vmpi: Parameters and datatypes
 * ***************************************************************************
 */


/** @brief Class Vmpi: Definition */

typedef struct Vmpi {
    int  mpi_rank;     /**< my process ID                                   */
    int  mpi_size;     /**< number of processess in this execution           */
} Vmpi;

/*
 * ***************************************************************************
 * Class Vmpi: Inlineable methods (vmpi.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_MALOC)
#else /* if defined(VINLINE_MALOC) */
#endif /* if !defined(VINLINE_MALOC) */

/** @brief Class Vmpi: Non-inlineable method (vmpi.c) */
VEXTERNC int Vmpi_init(int *argc, char ***argv);
/** @brief Class Vmpi: Non-inlineable method (vmpi.c) */
VEXTERNC int Vmpi_finalize(void);

/** @brief Class Vmpi: Non-inlineable method (vmpi.c) */
VEXTERNC Vmpi* Vmpi_ctor(void);
/** @brief Class Vmpi: Non-inlineable method (vmpi.c) */
VEXTERNC void Vmpi_dtor(Vmpi **thee);

/** @brief Class Vmpi: Non-inlineable method (vmpi.c) */
VEXTERNC int Vmpi_rank(Vmpi *thee);
/** @brief Class Vmpi: Non-inlineable method (vmpi.c) */
VEXTERNC int Vmpi_size(Vmpi *thee);
/** @brief Class Vmpi: Non-inlineable method (vmpi.c) */
VEXTERNC int Vmpi_barr(Vmpi *thee);

/** @brief Class Vmpi: Non-inlineable method (vmpi.c) */
VEXTERNC int Vmpi_send(Vmpi *thee, int des, char *buf, int bufsize);
/** @brief Class Vmpi: Non-inlineable method (vmpi.c) */
VEXTERNC int Vmpi_recv(Vmpi *thee, int src, char *buf, int bufsize);

/** @brief Class Vmpi: Non-inlineable method (vmpi.c) */
VEXTERNC int Vmpi_bcast(Vmpi *thee, char *buf, int bufsize);
/** @brief Class Vmpi: Non-inlineable method (vmpi.c) */
VEXTERNC int Vmpi_reduce(Vmpi *thee, char *sbuf, char *rbuf, int bufsize);
/** @brief Class Vmpi: Non-inlineable method (vmpi.c) */
VEXTERNC int Vmpi_isend(Vmpi *thee, int des, char *buf, int bufsize);

#endif /* _VMPI_H_ */

