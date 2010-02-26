/** @defgroup Vpbe Vpbe class
*  @brief    The Poisson-Boltzmann master class
*    
*            Contains objects and parameters used in every PBE calculation,
*            regardless of method.
* 
*/

/**
*  @file       vpbe.h
 *  @ingroup    Vpbe
 *  @brief      Contains declarations for class Vpbe
 *  @version  $Id: vpbe.h 1033 2006-12-29 17:08:22Z sobolevnrm $
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

#ifndef _VPBE_H_
#define _VPBE_H_

/* Generic headers */
#include "maloc/maloc.h"
#include "apbs/vhal.h"

/* Specific headers */
#include "apbs/vunit.h"
#include "apbs/vatom.h"
#include "apbs/vacc.h"
#include "apbs/vclist.h"

/**
*  @ingroup Vpbe
 *  @author  Nathan Baker
 *  @brief   Contains public data members for Vpbe class/module
 */
struct sVpbe { 
	
	Vmem *vmem;         /**< Memory management object */
	
	Valist *alist;      /**< Atom (charge) list */
	Vclist *clist;      /**< Atom location cell list */
	Vacc *acc;          /**< Accessibility object */
	
	double T;           /**< Temperature (K) */
	double soluteDiel;  /**< Solute dielectric constant (unitless) */
	double solventDiel; /**< Solvent dielectric constant (unitless) */
	double solventRadius;
	/**< Solvent probe radius (angstroms) for accessibility;
	* determining defining volumes for the dielectric
		* coefficient */
	double bulkIonicStrength; /**< Bulk ionic strength (M) */
	double maxIonRadius;      /**< Max ion radius (A; used for calculating
		* accessiblity and defining volumes for ionic
		* strength coeffcients) */
	int    numIon;            /**< Total number of ion species */
	double ionConc[MAXION];   /**< Concentration (M) of each species */
	double ionRadii[MAXION];  /**< Ionic radius (A) of each species */
	double ionQ[MAXION];      /**< Charge (e) of each species */
	
	double xkappa;      /**< Debye-Huckel parameter (bulk) */
	double deblen;      /**< Debye length (bulk) */
	double zkappa2;     /**< Square of modified Debye-Huckel parameter (bulk) */
	double zmagic;      /**< Delta function scaling parameter */
	
	double soluteCenter[3]; /**< Center of solute molecule (A) */
	double soluteRadius; /**< Radius of solute molecule (A) */
	double soluteXlen;  /**< Solute length in x-direction */
	double soluteYlen;  /**< Solute length in y-direction */
	double soluteZlen;  /**< Solute length in z-direction */
	double soluteCharge; /**< Charge of solute molecule (e) */
	
	int paramFlag;      /**< Check to see if the parameters have been set */
	
};

/** 
*  @ingroup Vpbe
*  @brief   Declaration of the Vpbe class as the Vpbe structure
*/
typedef struct sVpbe Vpbe;

/* ///////////////////////////////////////////////////////////////////////////
   // Class Vpbe: Inlineable methods (vpbe.c)
   /////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VPBE)

/** @brief   Get atom list
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee Vpbe object
*  @return  Pointer to internal Valist object
*/
VEXTERNC Valist* Vpbe_getValist(Vpbe *thee);

/** @brief   Get accessibility oracle
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee Vpbe object
*  @return  Pointer to internal Vacc object
*/
VEXTERNC Vacc*   Vpbe_getVacc(Vpbe *thee);

/** @brief   Get bulk ionic strength
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee Vpbe object
*  @return  Bulk ionic strength (M)
*/
VEXTERNC double  Vpbe_getBulkIonicStrength(Vpbe *thee);

/** @brief   Get maximum radius of ion species
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee Vpbe object
*  @return  Maximum radius (A)
*/
VEXTERNC double  Vpbe_getMaxIonRadius(Vpbe *thee);

/** @brief   Get temperature
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee Vpbe object
*  @return  Temperature (K)
*/
VEXTERNC double  Vpbe_getTemperature(Vpbe *thee);           

/** @brief   Get solute dielectric constant
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee Vpbe object
*  @return  Solute dielectric constant
*/
VEXTERNC double  Vpbe_getSoluteDiel(Vpbe *thee); 

/** @brief   Get apolar coefficient
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee Vpbe object
*  @return  Apolar coefficent (kJ/mol/A^2)
*/
VEXTERNC double  Vpbe_getGamma(Vpbe *thee);

/** @brief   Get sphere radius which bounds biomolecule
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee Vpbe object
*  @return  Sphere radius which bounds biomolecule (A)
*/
VEXTERNC double  Vpbe_getSoluteRadius(Vpbe *thee);

/** @brief   Get length of solute in x dimension
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee Vpbe object
*  @return  Length of solute in x dimension (A)
*/
VEXTERNC double  Vpbe_getSoluteXlen(Vpbe *thee);

/** @brief   Get length of solute in y dimension
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee Vpbe object
*  @return  Length of solute in y dimension (A)
*/
VEXTERNC double  Vpbe_getSoluteYlen(Vpbe *thee);

/** @brief   Get length of solute in z dimension
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee Vpbe object
*  @return  Length of solute in z dimension (A)
*/
VEXTERNC double  Vpbe_getSoluteZlen(Vpbe *thee);

/** @brief   Get coordinates of solute center
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee Vpbe object
*  @return  Pointer to 3*double array with solute center coordinates (A)
*/
VEXTERNC double* Vpbe_getSoluteCenter(Vpbe *thee);

/** @brief   Get total solute charge
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee Vpbe object
*  @return  Total solute charge (e)
*/
VEXTERNC double  Vpbe_getSoluteCharge(Vpbe *thee);

/** @brief   Get solvent dielectric constant
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee Vpbe object
*  @return  Solvent dielectric constant
*/
VEXTERNC double  Vpbe_getSolventDiel(Vpbe *thee);

/** @brief   Get solvent molecule radius
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee Vpbe object
*  @return  Solvent molecule radius (A)
*/
VEXTERNC double  Vpbe_getSolventRadius(Vpbe *thee);

/** @brief   Get Debye-Huckel parameter 
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee Vpbe object
*  @return  Bulk Debye-Huckel parameter (&Aring;)
*/
VEXTERNC double  Vpbe_getXkappa(Vpbe *thee);

/** @brief   Get Debye-Huckel screening length
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee Vpbe object
*  @return  Debye-Huckel screening length (&Aring;)
*/
VEXTERNC double  Vpbe_getDeblen(Vpbe *thee);

/** @brief   Get modified squared Debye-Huckel parameter
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee Vpbe object
*  @return  Modified squared Debye-Huckel parameter (\f$\AA^{-2}\f$)
*/
VEXTERNC double  Vpbe_getZkappa2(Vpbe *thee);

/** @brief   Get charge scaling factor
*  @ingroup Vpbe
*  @author  Nathan Baker and Mike Holst
*  @param   thee Vpbe object
*  @return  Get factor for scaling charges (in e) to internal units
*/
VEXTERNC double  Vpbe_getZmagic(Vpbe *thee);

#else /* if defined(VINLINE_VPBE) */
#   define Vpbe_getValist(thee) ((thee)->alist)
#   define Vpbe_getVacc(thee) ((thee)->acc)
#   define Vpbe_getBulkIonicStrength(thee) ((thee)->bulkIonicStrength)
#   define Vpbe_getTemperature(thee) ((thee)->T)           
#   define Vpbe_getSoluteDiel(thee) ((thee)->soluteDiel) 
#   define Vpbe_getSoluteCenter(thee) ((thee)->soluteCenter)
#   define Vpbe_getSoluteRadius(thee) ((thee)->soluteRadius)
#   define Vpbe_getSoluteXlen(thee) ((thee)->soluteXlen)
#   define Vpbe_getSoluteYlen(thee) ((thee)->soluteYlen)
#   define Vpbe_getSoluteZlen(thee) ((thee)->soluteZlen)
#   define Vpbe_getSoluteCharge(thee) ((thee)->soluteCharge)
#   define Vpbe_getSolventDiel(thee) ((thee)->solventDiel)
#   define Vpbe_getSolventRadius(thee) ((thee)->solventRadius)
#   define Vpbe_getMaxIonRadius(thee) ((thee)->maxIonRadius)
#   define Vpbe_getXkappa(thee) ((thee)->xkappa)
#   define Vpbe_getDeblen(thee) ((thee)->deblen)
#   define Vpbe_getZkappa2(thee) ((thee)->zkappa2)
#   define Vpbe_getZmagic(thee) ((thee)->zmagic)
#endif /* if !defined(VINLINE_VPBE) */

/* ///////////////////////////////////////////////////////////////////////////
   // Class Vpbe: Non-Inlineable methods (vpbe.c)
   /////////////////////////////////////////////////////////////////////////// */

/** @brief   Construct Vpbe object
*  @ingroup Vpbe
*  @author  Nathan Baker and Mike Holst
*  @note   This is partially based on some of Mike Holst's PMG code.  Here
*           are a few of the original function comments:
*           kappa is defined as follows:
*           \f[ \kappa^2 = \frac{8 \pi N_A e_c^2 I_s}{1000 \epsilon_w k_B T}
	*           \f] 
*           where the units are esu*esu/erg/mol.  To obtain \f$\AA^{-2}\f$, we
*           multiply by \f$10^{-16}\f$.  Thus, in \f$\AA^{-2}\f$, where
*           \f$k_B\f$ and \f$e_c\f$ are in gaussian rather than mks units, the
*           proper value for kappa is:
*           \f[ \kappa^2 = \frac{8 \pi N_A e_c^2 I_s}{1000 \epsilon_w k_b T}
	*           \times 10^{-16} \f]
*           and the factor of \f$10^{-16}\f$ results from converting cm^2 to
*           angstroms^2, noting that the 1000 in the denominator has converted
*           m^3 to cm^3, since the ionic strength \f$I_s\f$ is assumed to have
*           been provided in moles per liter, which is moles per 1000 cm^3.
*  @param   alist  Atom list
*  @param   ionNum  Number of counterion species
*  @param   ionConc Array containing counterion species' concentrations (M)
*  @param   ionRadii Array containing counterion species' radii (A)
*  @param   ionQ Array containing counterion species' charges (e)
*  @param   T temperature (K)
*  @param   soluteDiel Solute dielectric constant
*  @param   solventDiel Solvent dielectric constant
*  @param   solventRadius Solvent radius
*  @param   focusFlag 1 if Focusing operation, 0 otherwise
*  @param   sdens Vacc sphere density
*  @return  Pointer to newly allocated Vpbe object
*/

VEXTERNC Vpbe*   Vpbe_ctor(Valist *alist, int ionNum, double *ionConc, 
						   double *ionRadii, double *ionQ, double T, 
						   double soluteDiel, double solventDiel,  
						   double solventRadius, int focusFlag, double sdens);

/** @brief   FORTRAN stub to construct Vpbe objct
*  @ingroup Vpbe
*  @author  Nathan Baker and Mike Holst
*  @note   This is partially based on some of Mike Holst's PMG code.  Here
*           are a few of the original function comments:
*           kappa is defined as follows:
*           \f[ \kappa^2 = \frac{8 \pi N_A e_c^2 I_s}{1000 eps_w k_B T} \f]
*           where the units are esu*esu/erg/mol.  To obtain \f$\AA^{-2}\f$, we
*           multiply by \f$10^{-16}\f$.
*           Thus, in \f$\AA^{-2}\f$, where \f$k_B\f$ and \f$e_c\f$ are in 
*           gaussian rather than mks units, the proper value for kappa is:
*           \f[ \kappa^2 = \frac{8 pi N_A e_c^2 I_s}{1000 eps_w k_b T} \times 
	*           10^{-16} \f]
*           and the factor of \f$10^{-16}\f$ results from converting cm^2 to 
*           angstroms^2, noting that the 1000 in the denominator has converted
*           m^3 to cm^3, since the ionic strength \f$I_s\f$ is assumed to have
*           been provided in moles per liter, which is moles per 1000 cm^3. 
*  @param   thee   Pointer to memory allocated for Vpbe object
*  @param   alist  Atom list
*  @param   ionNum  Number of counterion species
*  @param   ionConc Array containing counterion species' concentrations (M)
*  @param   ionRadii Array containing counterion species' radii (A)
*  @param   ionQ Array containing counterion species' charges (e)
*  @param   T temperature (K)
*  @param   soluteDiel Solute dielectric constant
*  @param   solventDiel Solvent dielectric constant
*  @param   solventRadius Solvent radius
*  @param   focusFlag 1 if Focusing operation, 0 otherwise
*  @bug     The focusing flag is currently not used!!
*  @param   sdens Vacc sphere density
*  @return  1 if successful, 0 otherwise
*/
VEXTERNC int    Vpbe_ctor2(Vpbe *thee, Valist *alist, int ionNum, 
						   double *ionConc, double *ionRadii, double *ionQ, 
						   double T, double soluteDiel, 
						   double solventDiel, double solventRadius, int focusFlag,
						   double sdens);

/** @brief   Get information about the counterion species present
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee   Pointer to Vpbe object
*  @param   nion   Set to the number of counterion species
*  @param   ionConc Array to store counterion species' concentrations (M)
*  @param   ionRadii Array to store counterion species' radii (A)
*  @param   ionQ Array to store counterion species' charges (e)
*  @return  Number of ions
*/
VEXTERNC int     Vpbe_getIons(Vpbe *thee, int *nion, double ionConc[MAXION],
							  double ionRadii[MAXION], double ionQ[MAXION]);

/** @brief  Object destructor
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee   Pointer to memory location of object to be destroyed
*/
VEXTERNC void    Vpbe_dtor(Vpbe **thee);

/** @brief   FORTRAN stub object destructor
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee   Pointer to object to be destroyed
*/
VEXTERNC void    Vpbe_dtor2(Vpbe *thee);

/** @brief  Calculate coulombic energy of set of charges
*         
*           Perform an inefficient double sum to calculate the Coulombic
*           energy of a set of charges in a homogeneous dielectric (with
																	*           permittivity equal to the protein interior) and zero ionic
*           strength.  Result is returned in units of k_B T.  The sum can be
*           restriction to charges present in simplices of specified color
*           (pcolor); if (color == -1) no restrictions are used.
*
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee Vpbe object
*  @return  Coulombic energy in units of \f$k_B T\f$.
*/
VEXTERNC double  Vpbe_getCoulombEnergy1(Vpbe *thee);

/** @brief   Return the memory used by this structure (and its contents)
*           in bytes
*  @ingroup Vpbe
*  @author  Nathan Baker
*  @param   thee  Vpbe object
*  @return  The memory used by this structure and its contents in bytes
*/
VEXTERNC unsigned long int Vpbe_memChk(Vpbe *thee);

#endif /* ifndef _VPBE_H_ */
