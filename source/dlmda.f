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
c     epdtexp       polarization lambda exponent for dual topology
c     dedl          total energy lambda derivative
c     devdl         van der Waals lambda derivative
c     demdl         multipole lambda derivative
c     depdl         polarization lambda derivative
c     dedl2         total energy second order lambda derivative
c     devdl2        van der Waals second order lambda derivative
c     demdl2        multipole second order lambda derivative
c     depdl2        polarization second order lambda derivative
c     plambda       state weighting value for polarization potentials
c     dpldlmda      derivative of plambda wrt main lambda
c     deldlmda      derivative of elambda wrt main lambda
c     dvldlmda      derivative of vlambda wrt main lambda
c     d2pldlmda2    derivative of plambda wrt main lambda
c     d2eldlmda2    derivative of elambda wrt main lambda
c     d2vldlmda2    derivative of vlambda wrt main lambda
c     dvirdl        total virial lambda derivative
c     devvirdl      van der Waals virial lambda derivative
c     demvirdl      multipole virial lambda derivative
c     depvirdl      polarization virial lambda derivative
c     dfsumdl       total force lambda derivative
c     dfvdl         van der Waals force lambda derivative
c     dfmdl         multipole force lambda derivative
c     dfpdl         polarization force lambda derivative
c     pchgorig      original pchg
c     pchg0orig     original pchg0
c     bdplorig      original bdpl
c     poleorig      original pole
c     pcoreorig     original pcore
c     pvalorig      original pval
c     pval0orig     original pval0
c     polarityorig  original polarity
c     bflxorig      original bflx
c     aflxorig      original aflx
c     abflxorig     original abflx
c     lcmp          cmp for lambda derivative
c     lfmp          fmp for lambda derivative
c     lcphi         cphi for lambda derivative
c     lfphi         fphi for lambda derivative
c     lqgrid        qgrid for lambda derivative
c     use_dlmda     logical flag governing use of lambda derivative
c     douindorig    original douind
c
c
      module dlmda
      implicit none
      integer epdtexp
      real*8 dedl
      real*8 devdl
      real*8 demdl
      real*8 depdl
      real*8 dedl2
      real*8 devdl2
      real*8 demdl2
      real*8 depdl2
      real*8 plambda
      real*8 dpldlmda
      real*8 deldlmda
      real*8 dvldlmda
      real*8 d2pldlmda2
      real*8 d2eldlmda2
      real*8 d2vldlmda2
      real*8 dvirdl(3,3)
      real*8 devvirdl(3,3)
      real*8 demvirdl(3,3)
      real*8 depvirdl(3,3)
      real*8, allocatable :: dfsumdl(:,:)
      real*8, allocatable :: dfvdl(:,:)
      real*8, allocatable :: dfmdl(:,:)
      real*8, allocatable :: dfpdl(:,:)
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
