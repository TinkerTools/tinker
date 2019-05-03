c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  module mplpot  --  multipole functional form details  ##
c     ##                                                        ##
c     ############################################################
c
c
c     m2scale      scale factor for 1-2 multipole energy interactions
c     m3scale      scale factor for 1-3 multipole energy interactions
c     m4scale      scale factor for 1-4 multipole energy interactions
c     m5scale      scale factor for 1-5 multipole energy interactions
c     use_chgpen   flag to use charge penetration damped potential
c     pentyp       type of penetration damping (NONE, GORDON1, GORDON2)
c
c
      module mplpot
      implicit none
      real*8 m2scale,m3scale
      real*8 m4scale,m5scale
      logical use_chgpen
      character*7 pentyp
      save
      end
