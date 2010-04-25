c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  polpot.i  --  specifics of polarization functional form  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     poleps    induced dipole convergence criterion (rms Debyes/atom)
c     polsor    induced dipole SOR convergence acceleration factor
c     p2scale   field 1-2 scale factor for energy evaluations
c     p3scale   field 1-3 scale factor for energy evaluations
c     p4scale   field 1-4 scale factor for energy evaluations
c     p5scale   field 1-5 scale factor for energy evaluations
c     d1scale   field intra-group scale factor for direct induced
c     d2scale   field 1-2 group scale factor for direct induced
c     d3scale   field 1-3 group scale factor for direct induced
c     d4scale   field 1-4 group scale factor for direct induced
c     u1scale   field intra-group scale factor for mutual induced
c     u2scale   field 1-2 group scale factor for mutual induced
c     u3scale   field 1-3 group scale factor for mutual induced
c     u4scale   field 1-4 group scale factor for mutual induced
c     poltyp    type of polarization potential (direct or mutual)
c
c
      real*8 poleps,polsor
      real*8 p2scale,p3scale
      real*8 p4scale,p5scale
      real*8 d1scale,d2scale
      real*8 d3scale,d4scale
      real*8 u1scale,u2scale
      real*8 u3scale,u4scale
      character*6 poltyp
      common /polpot/ poleps,polsor,p2scale,p3scale,p4scale,p5scale,
     &                d1scale,d2scale,d3scale,d4scale,u1scale,u2scale,
     &                u3scale,u4scale,poltyp
