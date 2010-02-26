/**
 *  @file       vio.h
 *  @brief      Class Vio: virtual <SDIO/FILE/BUFF/UNIX/INET> I/O layer.
 *  @version    $Id: vio.h,v 1.21 2008/02/05 00:11:33 fetk Exp $
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


#ifndef _VIO_H_
#define _VIO_H_

#include <maloc/maloc_base.h>

#include <maloc/vnm.h>

/*
 * ***************************************************************************
 * Class Vio: Parameters and datatypes
 * ***************************************************************************
 */

#define VPORTNUMBER 14916   /**< our portbase;  5000 < VPORTNUMBER < 49152 */
#define VIO_MAXBUF 10       /**< number of internal buffers (BUFF datatype) */

/** @brief Class Vio: Parameters and datatypes */

typedef enum VIOtype {
    VIO_NO_TYPE, 
    VIO_SDIO,
    VIO_BUFF,
    VIO_FILE,
    VIO_UNIX,
    VIO_INET
} VIOtype;

/** @brief Class Vio: Parameters and datatypes */

typedef enum VIOfrmt {
    VIO_NO_FRMT,
    VIO_XDR,
    VIO_ASC
} VIOfrmt;

/** @brief Class Vio: Parameters and datatypes */

typedef enum VIOrwkey {
    VIO_NO_RW,
    VIO_R,
    VIO_W
} VIOrwkey;

/** @brief Class Vio: Definition */
typedef struct Vio {

    VIOtype type;       /**< file (or device) type.                          
                         *   VIO_NO_TYPE = not initialized.
                         *   VIO_SDIO    = standard I/O.
                         *   VIO_FILE    = file I/O.                   
                         *   VIO_BUFF    = buffer I/O.                      
                         *   VIO_UNIX    = UNIX (domain) socket I/O.
                         *   VIO_INET    = INET (network) socket I/O.       */

    VIOfrmt frmt;       /**< data format.
                         *   VIO_NO_FRMT = not initialized.                
                         *   VIO_ASC     = ASCII (FILE,BUFF,UNIX,INET).
                         *   VIO_XDR     = BINARY (FILE,BUFF,UNIX,INET).     */

    VIOrwkey rwkey;     /**< r/w key.
                         *   VIO_NO_R = not initialized.
                         *   VIO_R    = read (FILE,BUFF,UNIX,INET)
                         *   VIO_W    = write (FILE,BUFF,UNIX,INET)         */

    char file[VMAX_ARGLEN];   /**< file or device name (FILE,BUFF,UNIX,INET)  */
    char lhost[VMAX_ARGLEN];  /**< local hostname (me) (UNIX,INET)            */
    char rhost[VMAX_ARGLEN];  /**< remote hostname (other guy) (UNIX,INET)    */

    int error;          /**< note if any error has occurred on this vio device*/
    int dirty;          /**< dirty read bit -- have we read file yet (FILE)   */

    FILE *fp;           /**< file pointer (SDIO,FILE)                         */
    int so;             /**< primary unix domain or inet socket (UNIX,INET)   */
    int soc;            /**< subsocket created for socket reading (UNIX,INET) */
    void *name;         /**< &sockaddr_un or &sockaddr_in (UNIX,INET)         */
    void *axdr;         /**< ASC/XDR structure pointer (ASC,XDR)              */

    char whiteChars[VMAX_ARGNUM]; /**< white space character set (ASC)        */
    char commChars[VMAX_ARGNUM];  /**< comment character set (ASC,XDR)        */

    char ioBuffer[VMAX_BUFSIZE];  /**< I/O buffer (ASC,XDR)                   */
    int ioBufferLen;              /**< I/O buffer length (ASC,XDR)            */

    char putBuffer[VMAX_BUFSIZE]; /**< final write buffer (ASC,XDR)           */
    int putBufferLen;             /**< final write buffer length (ASC,XDR)    */

    char *VIObuffer;    /**< (BUFF) */
    int VIObufferLen;   /**< (BUFF) */
    int VIObufferPtr;   /**< (BUFF) */

} Vio;

/*
 * ***************************************************************************
 * Class Vio: Inlineable methods (vio.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_MALOC)
#else /* if defined(VINLINE_MALOC) */
#endif /* if !defined(VINLINE_MALOC) */


/** @brief Class Vio: Non-Inlineable method (vio.c) */
VEXTERNC void Vio_start(void);
/** @brief Class Vio: Non-Inlineable method (vio.c) */
VEXTERNC void Vio_stop(void);

/** @brief Class Vio: Non-Inlineable method (vio.c) */
VEXTERNC Vio* Vio_ctor(const char *socktype, const char *datafrmt, 
    const char *hostname, const char *filename, const char *rwkey);
/** @brief Class Vio: Non-Inlineable method (vio.c) */
VEXTERNC int Vio_ctor2(Vio *thee, const char *socktype, const char *datafrmt, 
    const char *hostname, const char *filename, const char *rwkey);

/** @brief Class Vio: Non-Inlineable method (vio.c) */
VEXTERNC void Vio_dtor(Vio **thee);
/** @brief Class Vio: Non-Inlineable method (vio.c) */
VEXTERNC void Vio_dtor2(Vio *thee);

/** @brief Class Vio: Non-Inlineable method (vio.c) */
VEXTERNC Vio *Vio_socketOpen(char *key,
    const char *iodev, const char *iofmt,
    const char *iohost, const char *iofile);
/** @brief Class Vio: Non-Inlineable method (vio.c) */
VEXTERNC void Vio_socketClose(Vio **sock);

/** @brief Class Vio: Non-Inlineable method (vio.c) */
VEXTERNC void Vio_setWhiteChars(Vio *thee, char *whiteChars);
/** @brief Class Vio: Non-Inlineable method (vio.c) */
VEXTERNC void Vio_setCommChars(Vio *thee, char *commChars);

/** @brief Class Vio: Non-Inlineable method (vio.c) */
VEXTERNC int Vio_accept(Vio *thee, int nonblock);
/** @brief Class Vio: Non-Inlineable method (vio.c) */
VEXTERNC void Vio_acceptFree(Vio *thee);

/** @brief Class Vio: Non-Inlineable method (vio.c) */
VEXTERNC int Vio_connect(Vio *thee, int nonblock);
/** @brief Class Vio: Non-Inlineable method (vio.c) */
VEXTERNC void Vio_connectFree(Vio *thee);

/** @brief Class Vio: Non-Inlineable method (vio.c) */
VEXTERNC int Vio_scanf(Vio *thee, char *parms, ...);
/** @brief Class Vio: Non-Inlineable method (vio.c) */
VEXTERNC int Vio_printf(Vio *thee, char *parms, ...);

/** @brief Class Vio: Non-Inlineable method (vio.c) */
VEXTERNC int Vio_read(Vio *thee, char *buf, int bufsize);
/** @brief Class Vio: Non-Inlineable method (vio.c) */
VEXTERNC int Vio_write(Vio *thee, char *buf, int bufsize);

/** @brief Class Vio: Non-Inlineable method (vio.c) */
VEXTERNC int Vio_bufSize(Vio *thee);
/** @brief Class Vio: Non-Inlineable method (vio.c) */
VEXTERNC char* Vio_bufGive(Vio *thee);
/** @brief Class Vio: Non-Inlineable method (vio.c) */
VEXTERNC void Vio_bufTake(Vio *thee, char *buf, int bufsize);

#endif /* _VIO_H_ */

