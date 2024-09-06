c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module math  --  mathematical and geometrical constants  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     pi         numerical value of the geometric constant Pi
c     twopi      numerical value of two times Pi
c     rootpi     numerical value of the square root of Pi
c     radian     conversion factor from radians to degrees
c     elog       numerical value of the natural logarithm base
c     logten     numerical value of the natural log of ten
c     twosix     numerical value of the sixth root of two
c     root2      numerical value of the square root of two
c     root3      numerical value of the square root of three
c     third      numerical value of one-third (1/3)
c     third2     numerical value of two-thirds (2/3)
c
c
      module math
      implicit none
      real*8 pi,twopi,rootpi
      real*8 radian,elog,logten
      real*8 twosix,root2,root3
      real*8 third,third2
      parameter (pi=3.141592653589793238d0)
      parameter (twopi=6.283185307179586476d0)
      parameter (rootpi=1.772453850905516027d0)
      parameter (radian=57.29577951308232088d0)
      parameter (elog=2.718281828459045235d0)
      parameter (logten=2.302585092994045684d0)
      parameter (twosix=1.122462048309372981d0)
      parameter (root2=1.414213562373095049d0)
      parameter (root3=1.732050807568877294d0)
      parameter (third=0.333333333333333333d0)
      parameter (third2=0.666666666666666667d0)
      save
      end
