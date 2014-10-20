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
c     maxtyp          force field atom type definitions
c     maxclass        force field atom class definitions
c     maxvlst         neighbors in van der Waals pair list
c     maxelst         neighbors in electrostatics pair list
c     maxulst         neighbors in dipole preconditioner list
c     maxval          atoms directly bonded to an atom
c     maxref          stored reference molecular systems
c     maxgrp          user-defined groups of atoms
c     maxres          residues in the macromolecule
c     maxfix          geometric constraints and restraints
c
c
      module sizes
      implicit none
      integer maxatm,maxtyp
      integer maxclass,maxvlst
      integer maxelst,maxulst
      integer maxval,maxref
      integer maxgrp,maxres
      integer maxfix
      parameter (maxatm=1000000)
      parameter (maxtyp=5000)
      parameter (maxclass=1000)
      parameter (maxvlst=1800)
      parameter (maxelst=1200)
      parameter (maxulst=100)
      parameter (maxval=8)
      parameter (maxref=30)
      parameter (maxgrp=1000)
      parameter (maxres=10000)
      parameter (maxfix=maxatm)
      save
      end
