c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  urypot.i  --  specifics of Urey-Bradley functional form  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     cury       cubic coefficient in Urey-Bradley potential
c     qury       quartic coefficient in Urey-Bradley potential
c     ureyunit   convert Urey-Bradley energy to kcal/mole
c
c
      real*8 cury,qury
      real*8 ureyunit
      common /urypot/ cury,qury,ureyunit
