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
 * rcsid="$Id: vmem.h,v 1.9 2002/10/01 21:29:45 mholst Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     vmem.h    < vmem.c >
 *
 * Purpose:  Class Vmem: A safer, object-oriented, malloc/free object.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#ifndef _VMEM_H_
#define _VMEM_H_

#include <maloc/maloc_base.h>

/*
 * ***************************************************************************
 * Class Vmem: Parameters and datatypes
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * Class Vmem: Definition
 * ***************************************************************************
 */

typedef struct Vmem {

    char name[80];      /* name of class we are managing malloc areas for   */

    int mallocBytes;    /* total size of all current malloc areas of class  */
    int freeBytes;      /* total size of all freed malloc areas of class    */
    int highWater;      /* high-water malloc bytemark for this class        */
    int mallocAreas;    /* total number of individual malloc areas          */

} Vmem;

/*
 * ***************************************************************************
 * Class Vmem: Inlineable methods (vmem.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_MALOC)
#else /* if defined(VINLINE_MALOC) */
#endif /* if !defined(VINLINE_MALOC) */

/*
 * ***************************************************************************
 * Class Vmem: Non-Inlineable methods (vmem.c)
 * ***************************************************************************
 */

VEXTERNC int Vmem_bytesTotal(void);
VEXTERNC int Vmem_mallocBytesTotal(void);
VEXTERNC int Vmem_freeBytesTotal(void);
VEXTERNC int Vmem_highWaterTotal(void);
VEXTERNC int Vmem_mallocAreasTotal(void);
VEXTERNC void Vmem_printTotal(void);

VEXTERNC Vmem* Vmem_ctor(char *name);
VEXTERNC void Vmem_dtor(Vmem **thee);

VEXTERNC void *Vmem_malloc(Vmem *thee, int num, int size);
VEXTERNC void Vmem_free(Vmem *thee, int num, int size, void **ram);
VEXTERNC void *Vmem_realloc(Vmem *thee, int num, int size, void **ram,
    int newNum);

VEXTERNC int Vmem_bytes(Vmem *thee);
VEXTERNC int Vmem_mallocBytes(Vmem *thee);
VEXTERNC int Vmem_freeBytes(Vmem *thee);
VEXTERNC int Vmem_highWater(Vmem *thee);
VEXTERNC int Vmem_mallocAreas(Vmem *thee);
VEXTERNC void Vmem_print(Vmem *thee);

#endif /* _VMEM_H_ */

