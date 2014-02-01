c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2005  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  dma.i  --  distributed multipole analysis components  ##
c     ##                                                        ##
c     ############################################################
c
c
c     mp        atomic monopole charge values from DMA
c     dpx       atomic dipole moment x-component from DMA
c     dpy       atomic dipole moment y-component from DMA
c     dpz       atomic dipole moment z-component from DMA
c     q20       atomic Q20 quadrupole component from DMA (zz)
c     q21c      atomic Q21c quadrupole component from DMA (xz)
c     q21s      atomic Q21s quadrupole component from DMA (yz)
c     q22c      atomic Q22c quadrupole component from DMA (xx-yy)
c     q22s      atomic Q22s quadrupole component from DMA (xy)
c
c
      real*8 mp,dpx,dpy,dpz
      real*8 q20,q21c,q21s
      real*8 q22c,q22s
      common /dma/ mp(maxatm),dpx(maxatm),dpy(maxatm),dpz(maxatm),
     &             q20(maxatm),q21c(maxatm),q21s(maxatm),q22c(maxatm),
     &             q22s(maxatm)
