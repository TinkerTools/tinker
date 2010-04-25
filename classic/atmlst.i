c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  atmlst.i  --  local geometry terms involving each atom  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     bndlist   list of the bond numbers involving each atom
c     anglist   list of the angle numbers centered on each atom
c
c
      integer bndlist,anglist
      common /atmlst/ bndlist(maxval,maxatm),
     &                anglist(maxval*(maxval-1)/2,maxatm)
