c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  katoms.i  --  forcefield parameters for the atom types  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     weight     average atomic mass of each atom type
c     atmcls     atom class number for each of the atom types
c     atmnum     atomic number for each of the atom types
c     ligand     number of atoms to be attached to each atom type
c     symbol     modified atomic symbol for each atom type
c     describe   string identifying each of the atom types
c
c
      integer atmcls,atmnum
      integer ligand
      real*8 weight
      character*3 symbol
      character*24 describe
      common /katoms/ weight(maxtyp),atmcls(maxtyp),atmnum(maxtyp),
     &                ligand(maxtyp),symbol(maxtyp),describe(maxtyp)
