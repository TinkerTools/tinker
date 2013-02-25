c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  sizes.i  --  parameter values to set array dimensions  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "sizes.i" sets values for critical array dimensions used
c     throughout the software; these parameters will fix the size
c     of the largest systems that can be handled; values too large
c     for the computer memory or swap space to accomodate will
c     result in poor performance or outright failure
c
c     parameter:      maximum allowed number of:
c
c     maxatm          atoms in the molecular system
c     maxval          atoms directly bonded to an atom
c     maxgrp          user-defined groups of atoms
c     maxref          stored reference molecular systems
c     maxtyp          force field atom type definitions
c     maxclass        force field atom class definitions
c     maxprm          lines in the parameter file
c     maxkey          lines in the keyword file
c     maxrot          bonds for torsional rotation
c     maxvar          optimization variables (vector storage)
c     maxopt          optimization variables (matrix storage)
c     maxhess         off-diagonal Hessian elements
c     maxlight        sites for method of lights neighbors
c     maxvlst         atom neighbors in van der Waals pair list
c     maxelst         atom neighbors in electrostatics pair list
c     maxfft          grid points in each FFT dimension
c     maxfix          geometric constraints and restraints
c     maxvib          vibrational frequencies
c     maxgeo          distance geometry points
c     maxcell         unit cells in replicated crystal
c     maxring         3-, 4-, or 5-membered rings
c     maxbio          biopolymer atom definitions
c     maxres          residues in the macromolecule
c     maxele          elements in periodic table
c     maxamino        amino acid residue types
c     maxnuc          nucleic acid residue types
c     maxbnd          covalent bonds in molecular system
c     maxang          bond angles in molecular system
c     maxtors         torsional angles in molecular system
c     maxbitor        bitorsions in molecular system
c
c
      integer maxatm,maxval,maxgrp
      integer maxref,maxtyp,maxclass
      integer maxprm,maxkey,maxrot
      integer maxvar,maxopt,maxhess
      integer maxlight,maxvlst,maxelst
      integer maxfft,maxfix,maxvib
      integer maxgeo,maxcell,maxring
      integer maxbio,maxres,maxele
      integer maxamino,maxnuc,maxbnd
      integer maxang,maxtors,maxbitor
      parameter (maxatm=1000000)
      parameter (maxval=8)
      parameter (maxgrp=1000)
      parameter (maxref=10)
      parameter (maxtyp=5000)
      parameter (maxclass=1000)
      parameter (maxprm=25000)
      parameter (maxkey=25000)
      parameter (maxrot=1000)
      parameter (maxvar=3*maxatm)
      parameter (maxopt=1000)
      parameter (maxhess=1000000)
      parameter (maxlight=8*maxatm)
      parameter (maxvlst=1800)
      parameter (maxelst=1200)
      parameter (maxfft=250)
      parameter (maxfix=maxatm)
      parameter (maxvib=1000)
      parameter (maxgeo=2500)
      parameter (maxcell=10000)
      parameter (maxring=10000)
      parameter (maxbio=10000)
      parameter (maxres=10000)
      parameter (maxele=112)
      parameter (maxamino=38)
      parameter (maxnuc=12)
      parameter (maxbnd=2*maxatm)
      parameter (maxang=4*maxatm)
      parameter (maxtors=6*maxatm)
      parameter (maxbitor=8*maxatm)
