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
c     elmdaexp      exponent for electrostatic exponential mapping
c     elmdainvn     inverse-power exponent for electrostatic mapping
c     emdtexp       multipole exponent for dual topology interpolation
c     epdtexp       polarization lambda exponent for dual topology
c     evdtexp       van der Waals exponent for dual topo interpolation
c     plmdaexp      exponent for polarization exponential mapping
c     plmdainvn     inverse-power exponent for polarization mapping
c     vlmdaexp      exponent for van der Waals exponential mapping
c     vlmdainvn     inverse-power exponent for van der Waals mapping
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
c     elmdainveps   shift for electrostatic inverse-power mapping
c     plmdainveps   shift for polarization inverse-power mapping
c     qntelmda0     sublambda lower bound for electrostatics
c     qntelmda1     sublambda upper bound for electrostatics
c     qntplmda0     sublambda lower bound for polarization
c     qntplmda1     sublambda upper bound for polarization
c     qntvlmda0     sublambda lower bound for van der Waals
c     qntvlmda1     sublambda upper bound for van der Waals
c     relstg1lmda0  main lambda where ligand 1 electrostatics start
c     relstg1lmda1  main lambda where ligand 1 electrostatics end
c     relstg2lmda0  main lambda where ligand 2 electrostatics start
c     relstg2lmda1  main lambda where ligand 2 electrostatics end
c     vlmdainveps   shift for van der Waals inverse-power mapping
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
c     relstagemix   flag that the active leg mixes two endpoint states
c     use_dlmda     logical flag governing use of lambda derivative
c     use_ele4f     flag to compute electrostatic final endpoint
c     use_ele4i     flag to compute electrostatic initial endpoint
c     use_emdt      flag governing use of multipole dual topology
c     use_evdt      flag governing use of van der Waals dual topology
c     use_epdt      flag governing use of polarization dual topology
c     use_meta      flag to use metadynamics
c     use_metadyn   flag to propagate metadynamics lambda particle
c     use_ost       flag to use orthogonal space tempering
c     use_ostdyn    flag to propagate lambda particle
c     use_plmda     flag governing rescale to a decoupled plambda
c     use_pol4f     flag to compute polarization final endpoint
c     use_pol4i     flag to compute polarization initial endpoint
c     use_relstage  flag to use staged relative free energy schedule
c     use_ti        flag to use thermodynamic integration
c     use_vdw4f     flag to compute van der Waals final endpoint
c     use_vdw4i     flag to compute van der Waals initial endpoint
c     douindorig    original douind
c     elmdamap      mapping type from main to electrostatic lambda
c     plmdamap      mapping type from main to polarization lambda
c     vlmdamap      mapping type from main to van der Waals lambda
c     lmdaengymode  free energy being computed, ABS or REL
c     lmdasampmode  method sampling the main lambda, META, OST or TI
c     relstage      active leg of the staged schedule
c
c
      module dlmda
      implicit none
      integer elmdaexp
      integer elmdainvn
      integer emdtexp
      integer epdtexp
      integer evdtexp
      integer plmdaexp
      integer plmdainvn
      integer vlmdaexp
      integer vlmdainvn
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
      real*8 elmdainveps
      real*8 plmdainveps
      real*8 qntelmda0
      real*8 qntelmda1
      real*8 qntplmda0
      real*8 qntplmda1
      real*8 qntvlmda0
      real*8 qntvlmda1
      real*8 relstg1lmda0
      real*8 relstg1lmda1
      real*8 relstg2lmda0
      real*8 relstg2lmda1
      real*8 vlmdainveps
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
      logical relstagemix
      logical use_dlmda
      logical use_ele4f
      logical use_ele4i
      logical use_emdt
      logical use_epdt
      logical use_evdt
      logical use_meta
      logical use_metadyn
      logical use_ost
      logical use_ostdyn
      logical use_plmda
      logical use_pol4f
      logical use_pol4i
      logical use_relstage
      logical use_ti
      logical use_vdw4f
      logical use_vdw4i
      logical, allocatable :: douindorig(:)
      character*3 elmdamap
      character*3 plmdamap
      character*3 vlmdamap
      character*4 lmdaengymode
      character*4 lmdasampmode
      character*4 relstage
      save
      end
