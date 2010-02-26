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
 * rcsid="$Id: vset.h,v 1.10 2002/10/01 21:29:45 mholst Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     vset.h    < vset.c >
 *
 * Purpose:  Class Vset: a dynamic set object.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#ifndef _VSET_H_
#define _VSET_H_

#include <maloc/maloc_base.h>

#include <maloc/vnm.h>
#include <maloc/vmem.h>

/*
 * ***************************************************************************
 * Class Vset: Parameters and datatypes
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * Class Vset: Definition
 * ***************************************************************************
 */

typedef struct Vset {

    Vmem *vmem;         /* the memory manager                               */
    int  iMadeVmem;     /* did i make vmem or was it inherited              */

    int curT;           /* the current "T" object in our collection         */

    char nameT[80];     /* name of object we are managing                   */
    int sizeT;          /* size of the object in bytes                      */

    int numBlocks;      /* total number of allocated blocks                 */
    int numT;           /* the global "T" counter -- how many "T"s in list  */
    int prtT;           /* for i/o at appropriate block creation/deletion   */

    int maxObjects;     /* number of objects to manage (user specified)     */
    int blockPower;     /* power of 2 for blocksize (e.g., =10, or =16)     */
    int blockSize;      /* blocksize is 2^(blockPower) */
    int blockMax;       /* num blocks = blockMax=(maxObjects/blockSize)     */
    int blockModulo;    /* =blockSize-1; for determining which block fast   */

    char **table;       /* list of pointers to blocks of storage we manage  */

} Vset;

/*
 * ***************************************************************************
 * Class Vset: Inlineable methods (vset.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_MALOC)
    VEXTERNC int Vset_num(Vset *thee);
    VEXTERNC char *Vset_access(Vset *thee, int i);
    VEXTERNC char *Vset_create(Vset *thee);
    VEXTERNC char *Vset_first(Vset *thee);
    VEXTERNC char *Vset_last(Vset *thee);
    VEXTERNC char *Vset_next(Vset *thee);
    VEXTERNC char *Vset_prev(Vset *thee);
    VEXTERNC char *Vset_peekFirst(Vset *thee);
    VEXTERNC char *Vset_peekLast(Vset *thee);
    VEXTERNC void Vset_destroy(Vset *thee);
#else /* if defined(VINLINE_MALOC) */
#   define Vset_num(thee) ((thee)->numT)
#   define Vset_access(thee,i) ( \
        ((i >= 0) && (i < thee->numT)) \
        ? &((thee)->table[ (i)>>(thee)->blockPower                 ] \
                         [ (thee)->sizeT*((i)&(thee)->blockModulo) ]) \
        : VNULL \
    )
#   define Vset_create(thee) ( \
        (  ((((thee)->numT)>>(thee)->blockPower) >= (thee)->numBlocks) \
        || ((((thee)->numT+1)%(thee)->prtT) == 0) ) \
        ? (Vset_createLast((thee))) \
        : (++((thee)->numT), (Vset_access((thee),(thee)->numT-1))) \
    )
#   define Vset_first(thee) ( \
        (thee)->curT = 0, \
        Vset_access((thee), (thee)->curT) \
    )
#   define Vset_last(thee) ( \
        (thee)->curT = (thee)->numT-1, \
        Vset_access((thee), (thee)->curT) \
    )
#   define Vset_next(thee) ( \
        (thee)->curT++, \
        ((thee)->curT < (thee)->numT) \
        ? Vset_access((thee), (thee)->curT) \
        : VNULL \
    )
#   define Vset_prev(thee) ( \
        (thee)->curT--, \
        ((thee)->curT >= 0) \
        ? Vset_access((thee), (thee)->curT) \
        : VNULL \
    )
#   define Vset_peekFirst(thee) ( \
        Vset_access((thee), 0) \
    )
#   define Vset_peekLast(thee) ( \
        Vset_access((thee), (thee)->numT-1) \
    )
#   define Vset_destroy(thee) ( \
        ( ((((thee)->numT-1)>>(thee)->blockPower) < (thee)->numBlocks-1) \
          || ((thee)->numT == 1) || ((((thee)->numT)%(thee)->prtT) == 0) ) \
        ? (Vset_destroyLast((thee))) : (void)(((thee)->numT)--) \
    )
#endif /* if !defined(VINLINE_MALOC) */

/*
 * ***************************************************************************
 * Class Vset: Non-Inlineable methods (vset.c)
 * ***************************************************************************
 */

VEXTERNC Vset* Vset_ctor(Vmem *vmem,
    const char *tname, int tsize, int tmaxNum, int ioKey);
VEXTERNC void Vset_dtor(Vset **thee);

VEXTERNC char *Vset_createLast(Vset *thee);
VEXTERNC void Vset_destroyLast(Vset *thee);
VEXTERNC void Vset_initData(Vset *thee);
VEXTERNC void Vset_reset(Vset *thee);
VEXTERNC void Vset_check(Vset *thee,
    int *tnum, int *tsize, int *tVecUse, int *tVecMal, int *tVecOhd);

VEXTERNC void Vset_memChk(Vset *thee);

#endif /* _VSET_H_ */

