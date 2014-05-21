c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  module sizes  --  parameters to set array dimensions  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "sizes" sets values for critical array dimensions used
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
c     maxrot          bonds for torsional rotation
c     maxopt          optimization variables (matrix storage)
c     maxvlst         neighbors in van der Waals pair list
c     maxelst         neighbors in electrostatics pair list
c     maxulst         neighbors in dipole preconditioner list
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
c
c
      module sizes
      implicit none
      integer maxatm,maxval
      integer maxgrp,maxref
      integer maxtyp,maxclass
      integer maxrot,maxopt
      integer maxvlst,maxelst
      integer maxulst,maxfft
      integer maxfix,maxvib
      integer maxgeo,maxcell
      integer maxring,maxbio
      integer maxres,maxele
      integer maxamino,maxnuc
      parameter (maxatm=1000000)
      parameter (maxval=8)
      parameter (maxgrp=1000)
      parameter (maxref=10)
      parameter (maxtyp=5000)
      parameter (maxclass=1000)
      parameter (maxrot=1000)
      parameter (maxopt=1000)
      parameter (maxvlst=1800)
      parameter (maxelst=1200)
      parameter (maxulst=100)
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
      save
      end
