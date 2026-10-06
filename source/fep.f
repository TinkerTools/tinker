c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  module fep  --  free energy perturbation variables  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     fepintv    dynamics steps between samples of the trial energies
c     nfep       total number of energy values recorded so far
c     nfepsave   total number of energy values already written out
c     nfeptot    total number of energy values the run can record
c     fepstep    dynamics step at which each energy was recorded
c     fepene     potential energy at the trial lambda of each record
c     feplmda    main lambda of the simulation for each record
c     feptrial   trial main lambda at which each energy was found
c
c
      module fep
      implicit none
      integer fepintv
      integer nfep
      integer nfepsave
      integer nfeptot
      integer, allocatable :: fepstep(:)
      real*8, allocatable :: fepene(:)
      real*8, allocatable :: feplmda(:)
      real*8, allocatable :: feptrial(:)
      save
      end
