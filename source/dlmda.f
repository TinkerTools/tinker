c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses Chung, Pengyu Ren, Jay Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ######################################################
c     ##                                                  ##
c     ##  module dlmda  --  lambda derivative components  ##
c     ##                                                  ##
c     ######################################################
c
c
c     epdtexp         polarization lambda exponent for dual topology
c     delmda          total energy lambda derivative
c     devlmda         van der Waals lambda derivative
c     demlmda         multipole lambda derivative
c     deplmda         polarization lambda derivative
c     delmda2         total energy second order lambda derivative
c     devlmda2        van der Waals second order lambda derivative
c     demlmda2        multipole second order lambda derivative
c     deplmda2        polarization second order lambda derivative
c     eplambda        state weighting value for polarization potentials
c     dldesum         total force lambda derivative
c     dldev           van der Waals force lambda derivative
c     dldem           multipole force lambda derivative
c     dldep           polarization force lambda derivative
c     pchgorig        original pchg
c     pchg0orig       original pchg0
c     bdplorig        original bdpl
c     poleorig        original pole
c     pcoreorig       original pcore
c     pvalorig        original pval
c     pval0orig       original pval0
c     polarityorig    original polarity
c     bflxorig        original bflx
c     aflxorig        original aflx
c     abflxorig       original abflx
c     lcmp            cmp for lambda derivative
c     lfmp            fmp for lambda derivative
c     lcphi           cphi for lambda derivative
c     lfphi           fphi for lambda derivative
c     lqgrid          qgrid for lambda derivative
c     use_dlmda       logical flag governing use of lambda derivative
c     douindorig      original douind
c
c
      module dlmda
      implicit none
      integer epdtexp
      real*8 delmda
      real*8 devlmda
      real*8 demlmda
      real*8 deplmda
      real*8 delmda2
      real*8 devlmda2
      real*8 demlmda2
      real*8 deplmda2
      real*8 eplambda
      real*8, allocatable :: dldesum(:,:)
      real*8, allocatable :: dldev(:,:)
      real*8, allocatable :: dldem(:,:)
      real*8, allocatable :: dldep(:,:)
      real*8, allocatable :: pchgorig(:)
      real*8, allocatable :: pchg0orig(:)
      real*8, allocatable :: bdplorig(:)
      real*8, allocatable :: poleorig(:,:)
      real*8, allocatable :: pcoreorig(:)
      real*8, allocatable :: pvalorig(:)
      real*8, allocatable :: pval0orig(:)
      real*8, allocatable :: polarityorig(:)
      real*8, allocatable :: bflxorig(:)
      real*8, allocatable :: aflxorig(:,:)
      real*8, allocatable :: abflxorig(:,:)
      real*8, allocatable :: lcmp(:,:)
      real*8, allocatable :: lfmp(:,:)
      real*8, allocatable :: lcphi(:,:)
      real*8, allocatable :: lfphi(:,:)
      real*8, allocatable :: lqgrid(:,:,:,:)
      logical use_dlmda
      logical, allocatable :: douindorig(:)
      save
      end
