c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
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
c     plambda         state weighting value for polarization potentials
c     dplambda        derivative of plambda wrt main lambda
c     delambda        derivative of elambda wrt main lambda
c     dvlambda        derivative of vlambda wrt main lambda
c     d2plambda       derivative of plambda wrt main lambda
c     d2elambda       derivative of elambda wrt main lambda
c     d2vlambda       derivative of vlambda wrt main lambda
c     dldvir          total virial lambda derivative
c     dldevvir        van der Waals virial lambda derivative
c     dldemvir        multipole virial lambda derivative
c     dldepvir        polarization virial lambda derivative
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
      real*8 plambda
      real*8 dplambda
      real*8 delambda
      real*8 dvlambda
      real*8 d2plambda
      real*8 d2elambda
      real*8 d2vlambda
      real*8 dldvir(3,3)
      real*8 dldevvir(3,3)
      real*8 dldemvir(3,3)
      real*8 dldepvir(3,3)
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
