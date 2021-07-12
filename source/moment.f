c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module moment  --  electric multipole moment components  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     netchg   net electric charge for the total system
c     netdpl   dipole moment magnitude for the total system
c     netqpl   diagonal quadrupole (Qxx, Qyy, Qzz) for system
c     xdpl     dipole vector x-component in the global frame
c     ydpl     dipole vector y-component in the global frame
c     zdpl     dipole vector z-component in the global frame
c     xxqpl    quadrupole tensor xx-component in global frame
c     xyqpl    quadrupole tensor xy-component in global frame
c     xzqpl    quadrupole tensor xz-component in global frame
c     yxqpl    quadrupole tensor yx-component in global frame
c     yyqpl    quadrupole tensor yy-component in global frame
c     yzqpl    quadrupole tensor yz-component in global frame
c     zxqpl    quadrupole tensor zx-component in global frame
c     zyqpl    quadrupole tensor zy-component in global frame
c     zzqpl    quadrupole tensor zz-component in global frame
c
c
      module moment
      implicit none
      real*8 netchg,netdpl
      real*8 netqpl(3)
      real*8 xdpl,ydpl,zdpl
      real*8 xxqpl,xyqpl,xzqpl
      real*8 yxqpl,yyqpl,yzqpl
      real*8 zxqpl,zyqpl,zzqpl
      save
      end
