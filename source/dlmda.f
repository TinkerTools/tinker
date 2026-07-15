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
c     emdtexp       multipole exponent for dual topology interpolation
c     epdtexp       polarization lambda exponent for dual topology
c     evdtexp       van der Waals exponent for dual topology interpolation
c     d2edl2        total energy second order lambda derivative
c     d2eldlmda2    second derivative of elambda wrt main lambda
c     d2emdl2       multipole second order lambda derivative
c     d2epdl2       polarization second order lambda derivative
c     d2evdl2       van der Waals second order lambda derivative
c     d2pldlmda2    second derivative of plambda wrt main lambda
c     d2vldlmda2    second derivative of vlambda wrt main lambda
c     dedl          total energy lambda derivative
c     deldlmda      derivative of elambda wrt main lambda
c     demdl         multipole lambda derivative
c     depdl         polarization lambda derivative
c     devdl         van der Waals lambda derivative
c     dpldlmda      derivative of plambda wrt main lambda
c     dvldlmda      derivative of vlambda wrt main lambda
c     plambda       state weighting value for polarization potentials
c     demvirdl      multipole virial lambda derivative
c     depvirdl      polarization virial lambda derivative
c     devvirdl      van der Waals virial lambda derivative
c     dvirdl        total virial lambda derivative
c     abflxorig     original abflx
c     aflxorig      original aflx
c     bdplorig      original bdpl
c     bflxorig      original bflx
c     dfmdl         multipole force lambda derivative
c     dfpdl         polarization force lambda derivative
c     dfsumdl       total force lambda derivative
c     dfvdl         van der Waals force lambda derivative
c     lcmp          cmp for lambda derivative
c     lcphi         cphi for lambda derivative
c     lfmp          fmp for lambda derivative
c     lfphi         fphi for lambda derivative
c     lqgrid        qgrid for lambda derivative
c     pchg0orig     original pchg0
c     pchgorig      original pchg
c     pcoreorig     original pcore
c     polarityorig  original polarity
c     poleorig      original pole
c     pval0orig     original pval0
c     pvalorig      original pval
c     use_dlmda     logical flag governing use of lambda derivative
c     use_emdt      flag governing use of multipole dual topology
c     use_evdt      flag governing use of van der Waals dual topology
c     use_epdt      flag governing use of polarization dual topology
c     douindorig    original douind
c
c
      module dlmda
      implicit none
      integer emdtexp
      integer epdtexp
      integer evdtexp
      real*8 d2edl2
      real*8 d2eldlmda2
      real*8 d2emdl2
      real*8 d2epdl2
      real*8 d2evdl2
      real*8 d2pldlmda2
      real*8 d2vldlmda2
      real*8 dedl
      real*8 deldlmda
      real*8 demdl
      real*8 depdl
      real*8 devdl
      real*8 dpldlmda
      real*8 dvldlmda
      real*8 plambda
      real*8 demvirdl(3,3)
      real*8 depvirdl(3,3)
      real*8 devvirdl(3,3)
      real*8 dvirdl(3,3)
      real*8, allocatable :: abflxorig(:,:)
      real*8, allocatable :: aflxorig(:,:)
      real*8, allocatable :: bdplorig(:)
      real*8, allocatable :: bflxorig(:)
      real*8, allocatable :: dfmdl(:,:)
      real*8, allocatable :: dfpdl(:,:)
      real*8, allocatable :: dfsumdl(:,:)
      real*8, allocatable :: dfvdl(:,:)
      real*8, allocatable :: lcmp(:,:)
      real*8, allocatable :: lcphi(:,:)
      real*8, allocatable :: lfmp(:,:)
      real*8, allocatable :: lfphi(:,:)
      real*8, allocatable :: lqgrid(:,:,:,:)
      real*8, allocatable :: pchg0orig(:)
      real*8, allocatable :: pchgorig(:)
      real*8, allocatable :: pcoreorig(:)
      real*8, allocatable :: polarityorig(:)
      real*8, allocatable :: poleorig(:,:)
      real*8, allocatable :: pval0orig(:)
      real*8, allocatable :: pvalorig(:)
      logical use_dlmda
      logical use_emdt
      logical use_epdt
      logical use_evdt
      logical, allocatable :: douindorig(:)
      save
      end
