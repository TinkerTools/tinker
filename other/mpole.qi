c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module mpole  --  atomic multipoles in current structure  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxpole   max components (monopole=1,dipole=4,quadrupole=13)
c
c     npole     total number of multipole sites in the system
c     ipole     number of the atom for each multipole site
c     polsiz    number of multipole components at each atom
c     pollist   multipole site for each atom (0=no multipole)
c     zaxis     number of the z-axis defining atom for each site
c     xaxis     number of the x-axis defining atom for each site
c     yaxis     number of the y-axis defining atom for each site
c     pole      multipole values for each site in the local frame
c     rpole     multipoles rotated to the global coordinate system
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     QI spherical harmonics stuff
c
c     spole     spherical multipole values in the local frame
c     srpole    spherical multipoles rotated to the global frame
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     polaxe    local axis type for each multipole site
c
c
      module mpole
      implicit none
      integer maxpole
      parameter (maxpole=13)
      integer npole
      integer, allocatable :: ipole(:)
      integer, allocatable :: polsiz(:)
      integer, allocatable :: pollist(:)
      integer, allocatable :: zaxis(:)
      integer, allocatable :: xaxis(:)
      integer, allocatable :: yaxis(:)
      real*8, allocatable :: pole(:,:)
      real*8, allocatable :: rpole(:,:)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     QI spherical harmonics stuff
c
c     real*8, allocatable :: spole(:,:)
c     real*8, allocatable :: srpole(:,:)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      character*8, allocatable :: polaxe(:)
      save
      end
