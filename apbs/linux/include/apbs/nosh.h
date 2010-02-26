/** @defgroup NOsh NOsh class
*  @brief    Class for parsing for fixed format input files
*/

/**
*  @file     nosh.h
 *  @ingroup  NOsh
 *  @brief    Contains declarations for class NOsh
 *  @version  $Id: nosh.h 1045 2007-01-02 16:16:42Z sdg0919 $
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

#ifndef _NOSH_H_
#define _NOSH_H_

/* Generic headers */
#include "maloc/maloc.h"
#include "apbs/vhal.h"

/* Headers specific to this file */
#include "apbs/pbeparm.h"
#include "apbs/mgparm.h"
#include "apbs/femparm.h"
#include "apbs/apolparm.h"
#include "apbs/valist.h"

/** @brief Maximum number of molecules in a run 
*  @ingroup NOsh */
#define NOSH_MAXMOL 20

/** @brief Maximum number of calculations in a run 
*  @ingroup NOsh */
#define NOSH_MAXCALC 20

/** @brief Maximum number of PRINT statements in a run 
*  @ingroup NOsh */
#define NOSH_MAXPRINT 20

/** @brief Maximum number of operations in a PRINT statement
*  @ingroup NOsh */
#define NOSH_MAXPOP 20

/**
* @brief  Molecule file format types
 * @ingroup NOsh
 */
enum eNOsh_MolFormat {
    NMF_PQR=0,  /**< PQR format */
    NMF_PDB=1,  /**< PDB format */
    NMF_XML=2   /**< XML format */
};

/**
* @brief  Declare NOsh_MolFormat type
 * @ingroup  NOsh
 */
typedef enum eNOsh_MolFormat NOsh_MolFormat;

/**
* @brief  NOsh calculation types
 * @ingroup NOsh
 */
enum eNOsh_CalcType {
    NCT_MG=0,  /**< Multigrid */
    NCT_FEM=1, /**< Finite element */
	NCT_APOL=2 /**< non-polar */
};

/**
* @brief  Declare NOsh_CalcType type
 * @ingroup  NOsh
 */
typedef enum eNOsh_CalcType NOsh_CalcType;

/**
* @brief  Parameter file format types
 * @ingroup NOsh
 */
enum eNOsh_ParmFormat {
    NPF_FLAT=0,  /**< Flat-file format */
    NPF_XML=1    /**< XML format */
};

/**
* @brief  Declare NOsh_ParmFormat type
 * @ingroup  NOsh
 */
typedef enum eNOsh_ParmFormat NOsh_ParmFormat;

/**
* @brief  NOsh print types
 * @ingroup  NOsh
 */
enum eNOsh_PrintType {
    NPT_ENERGY=0, /**< Energy (deprecated) */
    NPT_FORCE=1, /**< Force (deprecated) */
	NPT_ELECENERGY, /**< Elec Energy */
    NPT_ELECFORCE, /**< Elec Force */
	NPT_APOLENERGY, /**< Apol Energy */
    NPT_APOLFORCE /**< Apol Force */
};

/**
* @brief Declare NOsh_PrintType type
 * @ingroup NOsh
 */
typedef enum eNOsh_PrintType NOsh_PrintType;

/**
*  @ingroup NOsh
*  @author  Nathan Baker
*  @brief   Calculation class for use when parsing fixed format input files
*/
struct sNOsh_calc {
	MGparm *mgparm;         /**< Multigrid parameters */
	FEMparm *femparm;       /**< Finite element parameters */
	PBEparm *pbeparm;       /**< Generic PBE parameters */
	APOLparm *apolparm;		/**< Non-polar parameters */
	NOsh_CalcType calctype; /**< Calculation type */
};

/** @typedef NOsh_calc
*  @ingroup NOsh
*  @brief   Declaration of the NOsh_calc class as the NOsh_calc structure
*/
typedef struct sNOsh_calc NOsh_calc;

/**
*  @ingroup NOsh
*  @author  Nathan Baker
*  @brief   Class for parsing fixed format input files
*/
struct sNOsh {
	
	NOsh_calc *calc[NOSH_MAXCALC];  /**< The array of calculation objects 
		corresponding to actual calculations performed by the code.  Compare to 
		sNOsh::elec */
    int ncalc;  /**< The number of calculations in the calc array */
	
	NOsh_calc *elec[NOSH_MAXCALC];  /**< The array of calculation objects
		corresponding to ELEC statements read in the input file.  Compare to 
		sNOsh::calc */
    int nelec;  /**< The number of elec statements in the input file and in the 
		elec array */
	
	NOsh_calc *apol[NOSH_MAXCALC];  /**< The array of calculation objects
		corresponding to APOLAR statements read in the input file.  Compare to 
		sNOsh::calc */
    int napol;  /**< The number of apolar statements in the input file and in the 
		apolar array */
	
    int ispara;  /**< 1 => is a parallel calculation, 0 => is not */
    int proc_rank;  /**< Processor rank in parallel calculation */
    int proc_size;  /**< Number of processors in parallel calculation */
    int bogus;  /**< A flag which tells routines using NOsh that this particular 
		NOsh is broken -- useful for parallel focusing calculations where the 
		user gave us too many processors (1 => ignore this NOsh; 0 => this NOsh 
										  is OK) */
    int elec2calc[NOSH_MAXCALC];  /**< A mapping between ELEC statements which 
		appear in the input file and calc objects stored above.  Since we allow 
		both normal and focused  multigrid, there isn't a 1-to-1 correspondence 
		between ELEC statements and actual calcualtions.  This can really 
		confuse operations which work on specific calculations further down the
		road (like PRINT).  Therefore this array is the initial point of entry 
		for any calculation-specific operation.  It points to a specific entry 
		in the calc array. */
	int apol2calc[NOSH_MAXCALC];  /**< (see elec2calc) */
	
    int nmol;  /**< Number of molecules */
    char molpath[NOSH_MAXMOL][VMAX_ARGLEN];   /**< Paths to mol files */
    NOsh_MolFormat molfmt[NOSH_MAXMOL];  /**< Mol files formats */
	Valist *alist[NOSH_MAXMOL];  /**<  Molecules for calculation (can be used in
		setting mesh centers */
    int gotparm;  /**< Either have (1) or don't have (0) parm */
    char parmpath[VMAX_ARGLEN];   /**< Paths to parm file */
    NOsh_ParmFormat parmfmt;  /**< Parm file format */
    int ndiel;  /**< Number of dielectric maps */
    char dielXpath[NOSH_MAXMOL][VMAX_ARGLEN];  /**< Paths to x-shifted 
		dielectric map files */
    char dielYpath[NOSH_MAXMOL][VMAX_ARGLEN];  /**< Paths to y-shifted 
		dielectric map files */
    char dielZpath[NOSH_MAXMOL][VMAX_ARGLEN];  /**< Paths to z-shifted 
		dielectric map files */
    Vdata_Format dielfmt[NOSH_MAXMOL];  /**< Dielectric maps file formats */
    int nkappa;  /**< Number of kappa maps */
    char kappapath[NOSH_MAXMOL][VMAX_ARGLEN]; /**< Paths to kappa map files */
    Vdata_Format kappafmt[NOSH_MAXMOL];  /**< Kappa maps file formats */
    int ncharge;  /**< Number of charge maps */
    char chargepath[NOSH_MAXMOL][VMAX_ARGLEN];   /**< Paths to charge map files */
    Vdata_Format chargefmt[NOSH_MAXMOL];  /**< Charge maps fileformats */
    int nprint;  /**< How many print sections? */
    NOsh_PrintType printwhat[NOSH_MAXPRINT];  /**< What do we print:  \li 0 = 
		energy, \li 1 = force */
    int printnarg[NOSH_MAXPRINT];  /**< How many arguments in energy list */
    int printcalc[NOSH_MAXPRINT][NOSH_MAXPOP]; /**< ELEC id (see elec2calc) */
    int printop[NOSH_MAXPRINT][NOSH_MAXPOP];  /**< Operation id (0 = add, 1 = 
		subtract) */
    int parsed;  /**< Have we parsed an input file yet? */
    char elecname[NOSH_MAXCALC][VMAX_ARGLEN]; /**< Optional user-specified name 
		for ELEC statement */
	char apolname[NOSH_MAXCALC][VMAX_ARGLEN]; /**< Optional user-specified name 
		for APOLAR statement */
};

/** 
*  @ingroup NOsh
*  @brief   Declaration of the NOsh class as the NOsh structure
*/ 
typedef struct sNOsh NOsh;

/* ///////////////////////////////////////////////////////////////////////////
   // Class NOsh: Inlineable methods (mcsh.c)
   /////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_NOSH)
/** @brief    Returns path to specified molecule
*  @ingroup  NOsh
*  @author   Nathan Baker
*  @param    thee Pointer to NOsh object
*  @param    imol Molecule ID of interest
*  @returns  Path string
*/
VEXTERNC char* NOsh_getMolpath(NOsh *thee, int imol);

/** @brief    Returns path to specified x-shifted dielectric map
*  @ingroup  NOsh
*  @author   Nathan Baker
*  @param    thee Pointer to NOsh object
*  @param    imap Map ID of interest
*  @returns  Path string
*/
VEXTERNC char* NOsh_getDielXpath(NOsh *thee, int imap);

/** @brief    Returns path to specified y-shifted dielectric map
*  @ingroup  NOsh
*  @author   Nathan Baker
*  @param    thee Pointer to NOsh object
*  @param    imap Map ID of interest
*  @returns  Path string
*/
VEXTERNC char* NOsh_getDielYpath(NOsh *thee, int imap);

/** @brief    Returns path to specified z-shifted dielectric map
*  @ingroup  NOsh
*  @author   Nathan Baker
*  @param    thee Pointer to NOsh object
*  @param    imap Map ID of interest
*  @returns  Path string
*/
VEXTERNC char* NOsh_getDielZpath(NOsh *thee, int imap);

/** @brief    Returns path to specified kappa map
*  @ingroup  NOsh
*  @author   Nathan Baker
*  @param    thee Pointer to NOsh object
*  @param    imap Map ID of interest
*  @returns  Path string
*/
VEXTERNC char* NOsh_getKappapath(NOsh *thee, int imap);

/** @brief    Returns path to specified charge distribution map
*  @ingroup  NOsh
*  @author   Nathan Baker
*  @param    thee Pointer to NOsh object
*  @param    imap Map ID of interest
*  @returns  Path string
*/
VEXTERNC char* NOsh_getChargepath(NOsh *thee, int imap);

/** @brief    Returns specified calculation object
*  @ingroup  NOsh
*  @author   Nathan Baker
*  @param    thee Pointer to NOsh object
*  @param    icalc Calculation ID of interest
*  @returns  Pointer to specified calculation object
*/
VEXTERNC NOsh_calc* NOsh_getCalc(NOsh *thee, int icalc);

/** @brief    Returns format of specified dielectric map
*  @ingroup  NOsh
*  @author   Nathan Baker
*  @param    thee Pointer to NOsh object
*  @param    imap Calculation ID of interest
*  @returns  Format of dielectric map
*/
VEXTERNC int NOsh_getDielfmt(NOsh *thee, int imap);

/** @brief    Returns format of specified kappa map
*  @ingroup  NOsh
*  @author   Nathan Baker
*  @param    thee Pointer to NOsh object
*  @param    imap Calculation ID of interest
*  @returns  Format of kappa map
*/
VEXTERNC int NOsh_getKappafmt(NOsh *thee, int imap);

/** @brief    Returns format of specified charge map
*  @ingroup  NOsh
*  @author   Nathan Baker
*  @param    thee Pointer to NOsh object
*  @param    imap Calculation ID of interest
*  @returns  Format of charge map
*/
VEXTERNC int NOsh_getChargefmt(NOsh *thee, int imap);

#else

#   define NOsh_getMolpath(thee, imol) ((thee)->molpath[(imol)])
#   define NOsh_getDielXpath(thee, imol) ((thee)->dielXpath[(imol)])
#   define NOsh_getDielYpath(thee, imol) ((thee)->dielYpath[(imol)])
#   define NOsh_getDielZpath(thee, imol) ((thee)->dielZpath[(imol)])
#   define NOsh_getKappapath(thee, imol) ((thee)->kappapath[(imol)])
#   define NOsh_getChargepath(thee, imol) ((thee)->chargepath[(imol)])
#   define NOsh_getCalc(thee, icalc) ((thee)->calc[(icalc)])
#   define NOsh_getDielfmt(thee, imap) ((thee)->dielfmt[(imap)])
#   define NOsh_getKappafmt(thee, imap) ((thee)->kappafmt[(imap)])
#   define NOsh_getChargefmt(thee, imap) ((thee)->chargefmt[(imap)])

#endif


/* ///////////////////////////////////////////////////////////////////////////
   // Class NOsh: Non-inlineable methods (mcsh.c)
   /////////////////////////////////////////////////////////////////////////// */

/** @brief   Return an integer ID of the observable to print (@see printwhat)
*  @ingroup NOsh
*  @author  Nathan Baker
*  @param   thee NOsh object to use
*  @param   iprint ID of PRINT statement
*  @returns An integer ID of the observable to print (@see printwhat)
*/
VEXTERNC NOsh_PrintType NOsh_printWhat(NOsh *thee, int iprint);

/** @brief   Return an integer mapping of an ELEC statement to a calculation ID
*           (@see elec2calc)
*  @ingroup NOsh
*  @author  Nathan Baker
*  @param   thee NOsh object to use
*  @param   ielec ID of ELEC statement
*  @returns An integer mapping of an ELEC statement to a calculation ID
*           (@see elec2calc)
*/
VEXTERNC char* NOsh_elecname(NOsh *thee, int ielec);

/** @brief   Return the name of an elec statement        
*  @ingroup NOsh
*  @author  Todd Dolinsky
*  @param   thee NOsh object to use
*  @param   icalc ID of CALC statement
*  @returns The name (if present) of an ELEC statement
*/
VEXTERNC int NOsh_elec2calc(NOsh *thee, int icalc);

/** @brief   Return the name of an apol statement        
*  @ingroup NOsh
*  @author  David Gohara
*  @param   thee NOsh object to use
*  @param   icalc ID of CALC statement
*  @returns The name (if present) of an APOL statement
*/
VEXTERNC int NOsh_apol2calc(NOsh *thee, int icalc);

/** @brief   Return number of arguments to PRINT statement (@see printnarg)
*  @ingroup NOsh
*  @author  Nathan Baker
*  @param   thee NOsh object to use
*  @param   iprint ID of PRINT statement
*  @returns Number of arguments to PRINT statement (@see printnarg)
*/
VEXTERNC int NOsh_printNarg(NOsh *thee, int iprint);

/** @brief   Return integer ID for specified operation (@see printop)
*  @ingroup NOsh
*  @author  Nathan Baker
*  @param   thee NOsh object to use
*  @param   iprint ID of PRINT statement
*  @param   iarg ID of operation in PRINT statement
*  @returns Integer ID for specified operation (@see printop)
*/ 
VEXTERNC int NOsh_printOp(NOsh *thee, int iprint, int iarg);

/** @brief   Return calculation ID for specified PRINT statement 
*           (@see printcalc)
*  @ingroup NOsh
*  @author  Nathan Baker
*  @param   thee NOsh object to use
*  @param   iprint ID of PRINT statement
*  @param   iarg ID of operation in PRINT statement
*  @returns Calculation ID for specified PRINT statement 
*           (@see printcalc)
*/
VEXTERNC int NOsh_printCalc(NOsh *thee, int iprint, int iarg);

/** @brief   Construct NOsh
*  @ingroup NOsh
*  @author  Nathan Baker
*  @param   rank   Rank of current processor in parallel calculation (0 if not 
*				parallel)
*  @param   size   Number of processors in parallel calculation (1 if not
*				parallel)
*  @returns Newly allocated and initialized NOsh object
*/
VEXTERNC NOsh* NOsh_ctor(int rank, int size);

/**	@brief	Construct NOsh_calc
*	@ingroup	NOsh
*	@author	Nathan Baker
*	@param	calcType	Calculation type
*	@returns	Newly allocated and initialized NOsh object
*/
VEXTERNC NOsh_calc* NOsh_calc_ctor(
							  NOsh_CalcType calcType
							  );

/**	@brief	Copy NOsh_calc object into thee
*	@ingroup	NOsh
*	@author	Nathan Baker
*	@param  thee	Target object
*	@param	source	Source object
*/
VEXTERNC int NOsh_calc_copy(
							NOsh_calc *thee,
							NOsh_calc *source
							);

/** @brief   Object destructor
*  @ingroup NOsh
*  @author  Nathan Baker
*  @param   thee  Pointer to memory location of NOsh_calc object
*/
VEXTERNC void  NOsh_calc_dtor(NOsh_calc **thee);

/** @brief   FORTRAN stub to construct NOsh
*  @ingroup NOsh
*  @author  Nathan Baker
*  @param   thee  Space for NOsh objet
*  @param   rank   Rank of current processor in parallel calculation (0 if not
*				parallel)
*  @param   size   Number of processors in parallel calculation (1 if not
*				parallel)
*  @returns 1 if successful, 0 otherwise
*/
VEXTERNC int   NOsh_ctor2(NOsh *thee, int rank, int size);

/** @brief   Object destructor
*  @ingroup NOsh
*  @author  Nathan Baker
*  @param   thee  Pointer to memory location of NOsh object
*/
VEXTERNC void  NOsh_dtor(NOsh **thee);

/** @brief   FORTRAN stub for object destructor
*  @ingroup NOsh
*  @author  Nathan Baker
*  @param   thee  Pointer to NOsh object
*/
VEXTERNC void  NOsh_dtor2(NOsh *thee);

/**	@brief   Parse an input file from a socket
*	@ingroup	NOsh
*	@note	Should be called before NOsh_setupCalc
*	@author		Nathan Baker and Todd Dolinsky
*	@param	thee  Pointer to NOsh object
*	@param	sock  Stream of tokens to parse
*	@return  1 if successful, 0 otherwise
*/
VEXTERNC int   NOsh_parseInput(NOsh *thee, Vio *sock);

/**	@brief	Parse an input file only from a file
*	@note	Included for SWIG wrapper compatibility
*	@note	Should be called before NOsh_setupCalc
*	@ingroup NOsh
*	@author  Nathan Baker and Todd Dolinsky
*	@param   thee      Pointer to NOsh object
*	@param   filename  Name/path of readable file
*	@return  1 if successful, 0 otherwise
*/
VEXTERNC int   NOsh_parseInputFile(NOsh *thee, char *filename);

/**	@brief	Setup the series of electrostatics calculations
*	@note	Should be called after NOsh_parseInput*
*	@ingroup	NOsh
*	@author	Nathan Baker and Todd Dolinsky
*	@param	thee	Pointer to NOsh object
*	@param	alist	Array of pointers to Valist objects (molecules used to center
														 mesh);
*	@return	1 if successful, 0 otherwise
*/
VEXTERNC int NOsh_setupElecCalc(
							NOsh *thee, /**< NOsh object */
							Valist *alist[NOSH_MAXMOL] /**< Atom list for calculation */
							);

/**	@brief	Setup the series of non-polar calculations
*	@note	Should be called after NOsh_parseInput*
*	@ingroup	NOsh
*	@author	Nathan Baker and Todd Dolinsky
*	@param	thee	Pointer to NOsh object
*	@param	alist	Array of pointers to Valist objects (molecules used to center
														 mesh);
*	@return	1 if successful, 0 otherwise
*/
VEXTERNC int NOsh_setupApolCalc(
								NOsh *thee, /**< NOsh object */
								Valist *alist[NOSH_MAXMOL] /**< Atom list for calculation */
								);

#endif 

