/**
 *  @file       vset.h
 *  @brief      Class Vset: a dynamic set object.
 *  @version    $Id: vset.h,v 1.13 2008/02/05 00:11:34 fetk Exp $
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

/** @brief Class Vset: Definition */
typedef struct Vset {

    Vmem *vmem;         /**< the memory manager                             */
    int  iMadeVmem;     /**< did i make vmem or was it inherited            */

    int curT;           /**< the current "T" object in our collection       */

    char nameT[VMAX_ARGLEN]; /**< name of object we are managing            */
    int sizeT;          /**< size of the object in bytes                    */

    int numBlocks;      /**< total number of allocated blocks               */
    int numT;           /**< the global "T" counter -- how many "T"s in list*/
    int prtT;           /**< for i/o at appropriate block creation/deletion */

    int maxObjects;     /**< number of objects to manage (user specified)   */
    int blockPower;     /**< power of 2 for blocksize (e.g., =10, or =16)   */
    int blockSize;      /**< blocksize is 2^(blockPower) */
    int blockMax;       /**< num blocks = blockMax=(maxObjects/blockSize)   */
    int blockModulo;    /**< =blockSize-1; for determining which block fast */

    char **table;       /**< list of pointers to blocks of storage we manage*/

} Vset;


#if !defined(VINLINE_MALOC)
    /** @brief Class Vset: Inlineable method (vset.c) */
    VEXTERNC int Vset_num(Vset *thee);
    /** @brief Class Vset: Inlineable method (vset.c) */
    VEXTERNC char *Vset_access(Vset *thee, int i);
    /** @brief Class Vset: Inlineable method (vset.c) */
    VEXTERNC char *Vset_create(Vset *thee);
    /** @brief Class Vset: Inlineable method (vset.c) */
    VEXTERNC char *Vset_first(Vset *thee);
    /** @brief Class Vset: Inlineable method (vset.c) */
    VEXTERNC char *Vset_last(Vset *thee);
    /** @brief Class Vset: Inlineable method (vset.c) */
    VEXTERNC char *Vset_next(Vset *thee);
    /** @brief Class Vset: Inlineable method (vset.c) */
    VEXTERNC char *Vset_prev(Vset *thee);
    /** @brief Class Vset: Inlineable method (vset.c) */
    VEXTERNC char *Vset_peekFirst(Vset *thee);
    /** @brief Class Vset: Inlineable method (vset.c) */
    VEXTERNC char *Vset_peekLast(Vset *thee);
    /** @brief Class Vset: Inlineable method (vset.c) */
    VEXTERNC void Vset_destroy(Vset *thee);
#else /* if defined(VINLINE_MALOC) */
    /** @brief Class Vset: Inlineable method (vset.c) */
#   define Vset_num(thee) ((thee)->numT)
    /** @brief Class Vset: Inlineable method (vset.c) */
#   define Vset_access(thee,i) ( \
        ((i >= 0) && (i < thee->numT)) \
        ? &((thee)->table[ (i)>>(thee)->blockPower                 ] \
                         [ (thee)->sizeT*((i)&(thee)->blockModulo) ]) \
        : VNULL \
    )
    /** @brief Class Vset: Inlineable method (vset.c) */
#   define Vset_create(thee) ( \
        (  ((((thee)->numT)>>(thee)->blockPower) >= (thee)->numBlocks) \
        || ((((thee)->numT+1)%(thee)->prtT) == 0) ) \
        ? (Vset_createLast((thee))) \
        : (++((thee)->numT), (Vset_access((thee),(thee)->numT-1))) \
    )
    /** @brief Class Vset: Inlineable method (vset.c) */
#   define Vset_first(thee) ( \
        (thee)->curT = 0, \
        Vset_access((thee), (thee)->curT) \
    )
    /** @brief Class Vset: Inlineable method (vset.c) */
#   define Vset_last(thee) ( \
        (thee)->curT = (thee)->numT-1, \
        Vset_access((thee), (thee)->curT) \
    )
    /** @brief Class Vset: Inlineable method (vset.c) */
#   define Vset_next(thee) ( \
        (thee)->curT++, \
        ((thee)->curT < (thee)->numT) \
        ? Vset_access((thee), (thee)->curT) \
        : VNULL \
    )
    /** @brief Class Vset: Inlineable method (vset.c) */
#   define Vset_prev(thee) ( \
        (thee)->curT--, \
        ((thee)->curT >= 0) \
        ? Vset_access((thee), (thee)->curT) \
        : VNULL \
    )
    /** @brief Class Vset: Inlineable method (vset.c) */
#   define Vset_peekFirst(thee) ( \
        Vset_access((thee), 0) \
    )
    /** @brief Class Vset: Inlineable method (vset.c) */
#   define Vset_peekLast(thee) ( \
        Vset_access((thee), (thee)->numT-1) \
    )
    /** @brief Class Vset: Inlineable method (vset.c) */
#   define Vset_destroy(thee) ( \
        ( ((((thee)->numT-1)>>(thee)->blockPower) < (thee)->numBlocks-1) \
          || ((thee)->numT == 1) || ((((thee)->numT)%(thee)->prtT) == 0) ) \
        ? (Vset_destroyLast((thee))) : (void)(((thee)->numT)--) \
    )
#endif /* if !defined(VINLINE_MALOC) */

 /** @brief Class Vset: Non-Inlineable method (vset.c) */
VEXTERNC Vset* Vset_ctor(Vmem *vmem,
    const char *tname, int tsize, int tmaxNum, int ioKey);
 /** @brief Class Vset: Non-Inlineable method (vset.c) */
VEXTERNC void Vset_dtor(Vset **thee);

 /** @brief Class Vset: Non-Inlineable method (vset.c) */
VEXTERNC char *Vset_createLast(Vset *thee);
 /** @brief Class Vset: Non-Inlineable method (vset.c) */
VEXTERNC void Vset_destroyLast(Vset *thee);
 /** @brief Class Vset: Non-Inlineable method (vset.c) */
VEXTERNC void Vset_initData(Vset *thee);
 /** @brief Class Vset: Non-Inlineable method (vset.c) */
VEXTERNC void Vset_reset(Vset *thee);
 /** @brief Class Vset: Non-Inlineable method (vset.c) */
VEXTERNC void Vset_check(Vset *thee,
    int *tnum, int *tsize, int *tVecUse, int *tVecMal, int *tVecOhd);

 /** @brief Class Vset: Non-Inlineable method (vset.c) */
VEXTERNC void Vset_memChk(Vset *thee);

#endif /* _VSET_H_ */

