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
c     nuse   total number of active atoms in energy calculation
c     iuse   numbers of the atoms active in energy calculation
c     use    true if an atom is active, false if inactive
c
c
      integer nuse,iuse
      logical use
      common /usage/ nuse,iuse(maxatm),use(0:maxatm)
