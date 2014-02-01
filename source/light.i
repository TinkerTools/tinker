c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  light.i  --  indices for method of lights pair neighbors  ##
c     ##                                                            ##
c     ################################################################
c
c
c     nlight  total number of sites for method of lights calculation
c     kbx     low index of neighbors of each site in the x-sorted list
c     kby     low index of neighbors of each site in the y-sorted list
c     kbz     low index of neighbors of each site in the z-sorted list
c     kex     high index of neighbors of each site in the x-sorted list
c     key     high index of neighbors of each site in the y-sorted list
c     kez     high index of neighbors of each site in the z-sorted list
c     locx    pointer from x-sorted list into original interaction list
c     locy    pointer from y-sorted list into original interaction list
c     locz    pointer from z-sorted list into original interaction list
c     rgx     pointer from original interaction list into x-sorted list
c     rgy     pointer from original interaction list into y-sorted list
c     rgz     pointer from original interaction list into z-sorted list
c
c
      integer nlight,kbx,kby,kbz
      integer kex,key,kez
      integer locx,locy,locz
      integer rgx,rgy,rgz
      common /light/ nlight,kbx(maxatm),kby(maxatm),kbz(maxatm),
     &               kex(maxatm),key(maxatm),kez(maxatm),
     &               locx(maxlight),locy(maxlight),locz(maxlight),
     &               rgx(maxlight),rgy(maxlight),rgz(maxlight)
