c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module katoms  --  atom definition forcefield parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     atmcls     atom class number for each of the atom types
c     atmnum     atomic number for each of the atom types
c     ligand     number of atoms to be attached to each atom type
c     weight     average atomic mass of each atom type
c     symbol     modified atomic symbol for each atom type
c     describe   string identifying each of the atom types
c
c
      module katoms
      use sizes
      implicit none
      integer atmcls(maxtyp)
      integer atmnum(maxtyp)
      integer ligand(maxtyp)
      real*8 weight(maxtyp)
      character*3 symbol(maxtyp)
      character*24 describe(maxtyp)
      save
      end
