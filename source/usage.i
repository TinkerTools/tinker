c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  usage.i  --  atoms active during energy computation  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     iuse   numbers of the atoms active in energy calculation
c     use    true if an atom is active, false if inactive
c     nuse   total number of active atoms in energy calculation
c
c
      integer nuse
      integer, pointer :: iuse(:)
      logical, pointer :: use(:)
      common /usage/ iuse,use,nuse
