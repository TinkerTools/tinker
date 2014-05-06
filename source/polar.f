c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module polar  --  induced dipole moments & polarizability  ##        
c     ##                                                             ##
c     #################################################################
c
c
c     npolar    total number of polarizable sites in the system
c     polarity  dipole polarizability for each multipole site (Ang**3)
c     thole     Thole polarizability damping value for each site
c     pdamp     value of polarizability scale factor for each site
c     uind      induced dipole components at each multipole site
c     uinp      induced dipoles in field used for energy interactions
c     uinds     GK or PB induced dipoles at each multipole site
c     uinps     induced dipoles in field used for GK or PB energy
c
c
      module polar
      implicit none
      integer npolar
      real*8, allocatable :: polarity(:)
      real*8, allocatable :: thole(:)
      real*8, allocatable :: pdamp(:)
      real*8, allocatable :: uind(:,:)
      real*8, allocatable :: uinp(:,:)
      real*8, allocatable :: uinds(:,:)
      real*8, allocatable :: uinps(:,:)
      save
      end
