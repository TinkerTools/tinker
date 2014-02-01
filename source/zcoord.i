c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  zcoord.i  --  Z-matrix internal coordinate definitions  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     zbond   bond length used to define each Z-matrix atom
c     zang    bond angle used to define each Z-matrix atom
c     ztors   angle or torsion used to define Z-matrix atom
c     iz      defining atom numbers for each Z-matrix atom
c
c
      integer iz
      real*8 zbond,zang,ztors
      common /zcoord/ zbond(maxatm),zang(maxatm),ztors(maxatm),
     &                iz(4,maxatm)
