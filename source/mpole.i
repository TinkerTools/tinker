c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  mpole.i  --  multipole components for current structure  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     maxpole   max components (monopole=1,dipole=4,quadrupole=13)
c
c     pole      multipole values for each site in the local frame
c     rpole     multipoles rotated to the global coordinate system
c     npole     total number of multipole sites in the system
c     ipole     number of the atom for each multipole site
c     polsiz    number of multipole components at each atom
c     pollist   multipole site for each atom (0=no multipole)
c     zaxis     number of the z-axis defining atom for each site
c     xaxis     number of the x-axis defining atom for each site
c     yaxis     number of the y-axis defining atom for each site
c     polaxe    local axis type for each multipole site
c
c
      integer maxpole
      parameter (maxpole=13)
      integer npole,ipole
      integer polsiz,pollist
      integer zaxis,xaxis,yaxis
      real*8 pole,rpole
      character*8 polaxe
      common /mpole/ pole(maxpole,maxatm),rpole(maxpole,maxatm),npole,
     &               ipole(maxatm),polsiz(maxatm),pollist(maxatm),
     &               zaxis(maxatm),xaxis(maxatm),yaxis(maxatm),
     &               polaxe(maxatm)
