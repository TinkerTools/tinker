/**
 *  @file       vcom.h
 *  @brief      Class Vcom: virtual (currently just MPI) communications layer
 *  @version    $Id: vcom.h,v 1.28 2008/02/05 00:10:01 fetk Exp $
 *  @authors    Nathan Baker and Michael Holst
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


#ifndef _VCOM_H_
#define _VCOM_H_

#include <maloc/maloc_base.h>

#include <maloc/vsys.h>

/** @brief A base value for MPI tags */
#define VCOM_MPI_TAG 111

/*
 * ***************************************************************************
 * Class Vcom: Parameters and datatypes
 * ***************************************************************************
 */

/**
 *  @brief   Class Vcom: Definition
 */
typedef struct Vcom {

    int  mpi_rank;     /**< Local PE rank from MPI                            */
    int  mpi_size;     /**< Total number of PEs in this communicator from MPI */ 

    int  type;         /**< Communications type.                              
                        *   0 = not initialized.                              
                        *   1 = Message Passing Interface 1.1                 */

    int  error;        /**< note if any error has occurred on this vcom device*/

    void *core;        /**< Private MPI core                                  */

} Vcom;

/*
 * ***************************************************************************
 * Class Vcom: Inlineable methods (vcom.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_MALOC)
#else /* if defined(VINLINE_MALOC) */
#endif /* if !defined(VINLINE_MALOC) */


/** @brief Class Vcom: Non-Inlineable methods (vcom.c) */
VEXTERNC int Vcom_init(int *argc, char ***argv);
/** @brief Class Vcom: Non-Inlineable methods (vcom.c) */
VEXTERNC int Vcom_finalize(void);

/** @brief Class Vcom: Non-Inlineable methods (vcom.c) */
VEXTERNC Vcom* Vcom_ctor(int commtype);
/** @brief Class Vcom: Non-Inlineable methods (vcom.c) */
VEXTERNC int Vcom_ctor2(Vcom* thee, int commtype);
/** @brief Class Vcom: Non-Inlineable methods (vcom.c) */
VEXTERNC void Vcom_dtor(Vcom **thee);
/** @brief Class Vcom: Non-Inlineable methods (vcom.c) */
VEXTERNC void Vcom_dtor2(Vcom *thee);

/** @brief Class Vcom: Non-Inlineable methods (vcom.c) */
VEXTERNC int Vcom_send(Vcom *thee, int des, void *buf, int len, int type, 
  int block);
/** @brief Class Vcom: Non-Inlineable methods (vcom.c) */
VEXTERNC int Vcom_recv(Vcom *thee, int src, void *buf, int len, int type, 
  int block);
/** @brief Class Vcom: Non-Inlineable methods (vcom.c) */
VEXTERNC int Vcom_getCount(Vcom *thee, int src, int *length, int type);
/** @brief Class Vcom: Non-Inlineable methods (vcom.c) */
VEXTERNC int Vcom_reduce(Vcom *thee, void *sendbuf, void *recvbuf, int length, 
  int type, int op);
/** @brief Class Vcom: Non-Inlineable methods (vcom.c) */
VEXTERNC int Vcom_size(Vcom *thee);
/** @brief Class Vcom: Non-Inlineable methods (vcom.c) */
VEXTERNC int Vcom_resize(Vcom *thee, int newsize);
/** @brief Class Vcom: Non-Inlineable methods (vcom.c) */
VEXTERNC int Vcom_rank(Vcom *thee);
/** @brief Class Vcom: Non-Inlineable methods (vcom.c) */
VEXTERNC int Vcom_barr(Vcom *thee);

#endif /* _VCOM_H_ */

