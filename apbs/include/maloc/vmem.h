/**
 *  @file       vmem.h
 *  @brief      Class Vmem: A safer, object-oriented, malloc/free object.
 *  @version    $Id: vmem.h,v 1.13 2008/03/12 05:13:59 fetk Exp $
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


#ifndef _VMEM_H_
#define _VMEM_H_

#include <maloc/maloc_base.h>

/*
 * ***************************************************************************
 * Class Vmem: Parameters and datatypes
 * ***************************************************************************
 */


/** @brief Class Vmem: Definition */
typedef struct Vmem {

    char name[VMAX_ARGLEN]; /**< name of class we manage malloc areas for   */

    size_t mallocBytes; /**< total size of all current malloc areas of class*/
    size_t freeBytes;   /**< total size of all freed malloc areas of class  */
    size_t highWater;   /**< high-water malloc bytemark for this class      */
    size_t mallocAreas; /**< total number of individual malloc areas        */

} Vmem;

/*
 * ***************************************************************************
 * Class Vmem: Inlineable methods (vmem.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_MALOC)
#else /* if defined(VINLINE_MALOC) */
#endif /* if !defined(VINLINE_MALOC) */

/** @brief Class Vmem: Non-Inlineable method (vmem.c) */
VEXTERNC size_t Vmem_bytesTotal(void);
/** @brief Class Vmem: Non-Inlineable method (vmem.c) */
VEXTERNC size_t Vmem_mallocBytesTotal(void);
/** @brief Class Vmem: Non-Inlineable method (vmem.c) */
VEXTERNC size_t Vmem_freeBytesTotal(void);
/** @brief Class Vmem: Non-Inlineable method (vmem.c) */
VEXTERNC size_t Vmem_highWaterTotal(void);
/** @brief Class Vmem: Non-Inlineable method (vmem.c) */
VEXTERNC size_t Vmem_mallocAreasTotal(void);
/** @brief Class Vmem: Non-Inlineable method (vmem.c) */
VEXTERNC void Vmem_printTotal(void);

/** @brief Class Vmem: Non-Inlineable method (vmem.c) */
VEXTERNC Vmem* Vmem_ctor(char *name);
/** @brief Class Vmem: Non-Inlineable method (vmem.c) */
VEXTERNC void Vmem_dtor(Vmem **thee);

/** @brief Class Vmem: Non-Inlineable method (vmem.c) */
VEXTERNC void *Vmem_malloc(Vmem *thee, size_t num, size_t size);
/** @brief Class Vmem: Non-Inlineable method (vmem.c) */
VEXTERNC void Vmem_free(Vmem *thee, size_t num, size_t size, void **ram);
/** @brief Class Vmem: Non-Inlineable method (vmem.c) */
VEXTERNC void *Vmem_realloc(Vmem *thee, size_t num, size_t size, void **ram,
    size_t newNum);

/** @brief Class Vmem: Non-Inlineable method (vmem.c) */
VEXTERNC size_t Vmem_bytes(Vmem *thee);
/** @brief Class Vmem: Non-Inlineable method (vmem.c) */
VEXTERNC size_t Vmem_mallocBytes(Vmem *thee);
/** @brief Class Vmem: Non-Inlineable method (vmem.c) */
VEXTERNC size_t Vmem_freeBytes(Vmem *thee);
/** @brief Class Vmem: Non-Inlineable method (vmem.c) */
VEXTERNC size_t Vmem_highWater(Vmem *thee);
/** @brief Class Vmem: Non-Inlineable method (vmem.c) */
VEXTERNC size_t Vmem_mallocAreas(Vmem *thee);
/** @brief Class Vmem: Non-Inlineable method (vmem.c) */
VEXTERNC void Vmem_print(Vmem *thee);

#endif /* _VMEM_H_ */

