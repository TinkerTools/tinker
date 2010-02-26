/** @defgroup Vstring Vstring class
 *  @brief    Provides a collection of useful non-ANSI string functions
 */

/**
 *  @file     vstring.h
 *  @ingroup  Vstring
 *  @brief    Contains declarations for class Vstring
 *  @version  $Id: vstring.h 1033 2006-12-29 17:08:22Z sobolevnrm $
 *  @author   Nathan A. Baker
 *  
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 * Nathan A. Baker (baker@biochem.wustl.edu)
 * Dept. of Biochemistry and Molecular Biophysics
 * Center for Computational Biology
 * Washington University in St. Louis
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2002-2007.  Washington University in St. Louis.
 * All Rights Reserved.
 * Portions Copyright (c) 1999-2002.  The Regents of the University of
 * California.  
 * Portions Copyright (c) 1995.  Michael Holst.
 *
 * This file is part of APBS.
 *
 * APBS is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * APBS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with APBS; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 *
 * Linking APBS statically or dynamically with other modules is making a
 * combined work based on APBS. Thus, the terms and conditions of the GNU
 * General Public License cover the whole combination.
 * 
 * SPECIAL GPL EXCEPTION
 * In addition, as a special exception, the copyright holders of APBS
 * give you permission to combine the APBS program with free software
 * programs and libraries that are released under the GNU LGPL or with
 * code included in releases of ISIM, Ion Simulator Interface, PMV, PyMOL
 * SMOL, VMD, and Vision. Such combined software may be linked with APBS and 
 * redistributed together in original or modified form as mere aggregation
 * without requirement that the entire work be under the scope of the GNU 
 * General Public License. This special exception permission is also extended
 * to any software listed in the SPECIAL GPL EXCEPTION clauses by the PMG,
 * FEtk, MC, or MALOC libraries.
 * 
 * Note that people who make modified versions of APBS are not obligated
 * to grant this special exception for their modified versions; it is
 * their choice whether to do so. The GNU General Public License gives
 * permission to release a modified version without this exception; this
 * exception also makes it possible to release a modified version which
 * carries forward this exception.
 *
 * @endverbatim
 */

#ifndef _VSTRING_H_
#define _VSTRING_H_

#include "maloc/maloc.h"
#include "apbs/vhal.h"

/** @brief   Case-insensitive string comparison (BSD standard)
 *  @ingroup Vstring
 *  @author  Copyright (c) 1988-1993 The Regents of the University of
 *           California.  Copyright (c) 1995-1996 Sun Microsystems, Inc.
 *  @note    Copyright (c) 1988-1993 The Regents of the University of
 *           California.  Copyright (c) 1995-1996 Sun Microsystems, Inc.
 *  @param   s1   First string for comparison
 *  @param   s2   Second string for comparison
 *  @return  An integer less than, equal to, or greater than zero if s1 is
 *           found,  respectively,  to  be  less  than, to match, or be greater
 *           than s2. (Source:  Linux man pages)
 */
VEXTERNC int Vstring_strcasecmp(const char *s1, const char *s2);

/** @brief   A modified sscanf that examines the complete string
 *  @ingroup Vstring
 *  @author  Todd Dolinsky
 *  @param   tok   The string to examine
 *  @return  1 if the entire string is an integer, 0 if otherwise.
 */
VEXTERNC int Vstring_isdigit(const char *tok);

#endif    /* ifndef _VSTRING_H_ */
