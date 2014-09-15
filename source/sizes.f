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
c     throughout the software; these parameters fix the size of
c     the largest systems that can be handled
c
c     parameter:      maximum allowed number of:
c
c     maxatm          atoms in the molecular system
c     maxval          atoms directly bonded to an atom
c     maxgrp          user-defined groups of atoms
c     maxref          stored reference molecular systems
c     maxtyp          force field atom type definitions
c     maxclass        force field atom class definitions
c     maxvlst         neighbors in van der Waals pair list
c     maxelst         neighbors in electrostatics pair list
c     maxulst         neighbors in dipole preconditioner list
c     maxfix          geometric constraints and restraints
c     maxres          residues in the macromolecule
c
c
      module sizes
      implicit none
      integer maxatm,maxval
      integer maxgrp,maxref
      integer maxtyp,maxclass
      integer maxvlst,maxelst
      integer maxulst,maxfix
      integer maxres
      parameter (maxatm=1000000)
      parameter (maxval=8)
      parameter (maxgrp=1000)
      parameter (maxref=10)
      parameter (maxtyp=5000)
      parameter (maxclass=1000)
      parameter (maxvlst=1800)
      parameter (maxelst=1200)
      parameter (maxulst=100)
      parameter (maxfix=maxatm)
      parameter (maxres=10000)
      save
      end
