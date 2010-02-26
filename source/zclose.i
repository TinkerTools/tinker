c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  zclose.i  --  ring openings and closures for Z-matrix  ##
c     ##                                                         ##
c     #############################################################
c
c
c     nadd   number of added bonds between Z-matrix atoms
c     iadd   numbers of the atom pairs defining added bonds
c     ndel   number of bonds between Z-matrix bonds to delete
c     idel   numbers of the atom pairs defining deleted bonds
c
c
      integer nadd,iadd
      integer ndel,idel
      common /zclose/ nadd,iadd(2,maxatm),ndel,idel(2,maxatm)
