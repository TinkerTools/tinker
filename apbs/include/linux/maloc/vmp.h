/**
 *  @file       vmp.h
 *  @brief      Class Vmp: a Virtual MPI communication layer object.
 *  @version    $Id: vmp.h,v 1.12 2008/02/05 00:10:01 fetk Exp $
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

/** @brief Class Vmp: Definition */
typedef struct Vmp {
    int  mpi_rank;     /**< my process ID                                   */
    int  mpi_size;     /**< number of processess in this execution          */
} Vmp;

/*
 * ***************************************************************************
 * Class Vmp: Inlineable methods (vmp.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_MALOC)
#else /* if defined(VINLINE_MALOC) */
#endif /* if !defined(VINLINE_MALOC) */


/** @brief Class Vmp: Non-inlineable method (vmp.c) */
VEXTERNC int Vmp_init(int *argc, char ***argv);
/** @brief Class Vmp: Non-inlineable method (vmp.c) */
VEXTERNC int Vmp_finalize(void);

/** @brief Class Vmp: Non-inlineable method (vmp.c) */
VEXTERNC Vmp* Vmp_ctor(void);
/** @brief Class Vmp: Non-inlineable method (vmp.c) */
VEXTERNC void Vmp_dtor(Vmp **thee);

/** @brief Class Vmp: Non-inlineable method (vmp.c) */
VEXTERNC int Vmp_rank(Vmp *thee);
/** @brief Class Vmp: Non-inlineable method (vmp.c) */
VEXTERNC int Vmp_size(Vmp *thee);
/** @brief Class Vmp: Non-inlineable method (vmp.c) */
VEXTERNC int Vmp_barr(Vmp *thee);

/** @brief Class Vmp: Non-inlineable method (vmp.c) */
VEXTERNC int Vmp_send(Vmp *thee, int des, char *buf, int bufsize);
/** @brief Class Vmp: Non-inlineable method (vmp.c) */
VEXTERNC int Vmp_recv(Vmp *thee, int src, char *buf, int bufsize);

#endif /* _VMP_H_ */

