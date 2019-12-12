c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2010  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module vibs  --  iterative vibrational analysis components  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     rho        trial vectors for iterative vibrational analysis
c     rhok       alternate vectors for iterative vibrational analysis
c     rwork      temporary work array for eigenvector transformation
c
c
      module vibs
      implicit none
      real*8, allocatable :: rho(:,:)
      real*8, allocatable :: rhok(:,:)
      real*8, allocatable :: rwork(:,:)
      save
      end
