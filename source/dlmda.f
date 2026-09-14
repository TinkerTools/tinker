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
c     elmdaapmn     power exponent for electrostatic asymmetric map
c     elmdaexp      exponent for electrostatic exponential mapping
c     elmdainvn     inverse-power exponent for electrostatic mapping
c     emdtexp       multipole exponent for dual topology interpolation
c     epdtexp       polarization lambda exponent for dual topology
c     erelst0       multipole coupling state at the lower endpoint
c     erelst1       multipole coupling state at the upper endpoint
c     evdtexp       van der Waals exponent for dual topo interpolation
c     lmdaintv      steps in each adaptive bias sample interval
c     lmdanpa       steps propagating the lambda particle
c     lmdanpb       steps equilibrating at the frozen lambda
c     lmdanpc       steps averaged at the frozen lambda
c     lmdastep      dynamics step count of the adaptive lambda bias
c     nlmda         number of lambda bins
c     nlmdahist     number of saved lambda bias history entries
c     nlmdasave     history entries already written to the file
c     nrelsub       number of parameter-zeroed subsystems, always five
c     plmdaapmn     power exponent for polarization asymmetric map
c     plmdaexp      exponent for polarization exponential mapping
c     plmdainvn     inverse-power exponent for polarization mapping
c     prelst0       polarization coupling state at the lower endpoint
c     prelst1       polarization coupling state at the upper endpoint
c     rellig1       coupling state id with ligand 1 bound to environment
c     rellig2       coupling state id with ligand 2 bound to environment
c     relnone       coupling state id with neither ligand bound
c     sizelmdahist  allocation size of the lambda bias history
c     vlmdaapmn     power exponent for van der Waals asymmetric map
c     vlmdaexp      exponent for van der Waals exponential mapping
c     vlmdainvn     inverse-power exponent for van der Waals mapping
c     vrelst0       van der Waals coupling state at the lower endpoint
c     vrelst1       van der Waals coupling state at the upper endpoint
c     lmdaihist     step at which each history entry was saved
c     d2edl2        total energy second order lambda derivative
c     d2eldlmda2    second derivative of elambda wrt main lambda
c     d2emdl2       multipole second order lambda derivative
c     d2epdl2       polarization second order lambda derivative
c     d2evdl2       van der Waals second order lambda derivative
c     d2pldlmda2    second derivative of plambda wrt main lambda
c     d2vldlmda2    second derivative of vlambda wrt main lambda
c     dedl          total unbiased energy lambda derivative
c     dedlavg       interval average of dU/dlambda
c     dedlstd       interval deviation of dU/dlambda
c     deffdl        effective lambda derivative for propagation
c     deldlmda      derivative of elambda wrt main lambda
c     demdl         multipole lambda derivative
c     depdl         polarization lambda derivative
c     devdl         van der Waals lambda derivative
c     dpldlmda      derivative of plambda wrt main lambda
c     dvldlmda      derivative of vlambda wrt main lambda
c     elmdaapmrho   endpoint slope ratio for electrostatic asym map
c     elmdainveps   shift for electrostatic inverse-power mapping
c     lmdaavg       interval average of the main lambda
c     lmdaddgdl     current dDeltaG/dlambda of the lambda bias
c     lmdadeltag    current free energy estimate of the lambda bias
c     lmdadfdl      saved bias free energy derivative for dynamics
c     lmdadt        time step of the theta lambda coordinate
c     lmdafric      friction coefficient of theta lambda coordinate
c     lmdamass      fictitious mass of the theta lambda coordinate
c     lmdaparatio   interval fraction propagating the lambda particle
c     lmdapbratio   interval fraction equilibrating at fixed lambda
c     lmdapcratio   interval fraction averaging at fixed lambda
c     lmdastd       interval deviation of the main lambda
c     lmdatheta     theta coordinate used to propagate lambda
c     lmdavbias     saved lambda bias energy shift
c     lmdavtheta    velocity of the theta lambda coordinate
c     plmdaapmrho   endpoint slope ratio for polarization asym map
c     plmdainveps   shift for polarization inverse-power mapping
c     qntelmda0     sublambda lower bound for electrostatics
c     qntelmda1     sublambda upper bound for electrostatics
c     qntplmda0     sublambda lower bound for polarization
c     qntplmda1     sublambda upper bound for polarization
c     qntvlmda0     sublambda lower bound for van der Waals
c     qntvlmda1     sublambda upper bound for van der Waals
c     vlmdaapmrho   endpoint slope ratio for van der Waals asym map
c     vlmdainveps   shift for van der Waals inverse-power mapping
c     wlmda         width of lambda bins
c     wlmda2        half width of lambda bins
c     demvirdl      multipole virial lambda derivative
c     depvirdl      polarization virial lambda derivative
c     devvirdl      van der Waals virial lambda derivative
c     dvirdl        total virial lambda derivative
c     bdplorig      original bdpl
c     bflxorig      original bflx
c     lmdafhist     dU/dlambda value of each history entry
c     lmdaflist     dU/dlambda values saved within an interval
c     lmdafmean     mean force of each lambda bin
c     lmdafsum      weighted dU/dlambda sum of each lambda bin
c     lmdafwt       total weight of each lambda bin
c     lmdalhist     lambda value of each history entry
c     lmdallist     lambda values saved within an interval
c     pchg0orig     original pchg0
c     pchgorig      original pchg
c     pcoreorig     original pcore
c     polarityorig  original polarity
c     pval0orig     original pval0
c     pvalorig      original pval
c     abflxorig     original abflx
c     aflxorig      original aflx
c     dfmdl         multipole force lambda derivative
c     dfpdl         polarization force lambda derivative
c     dfsumdl       total force lambda derivative
c     dfvdl         van der Waals force lambda derivative
c     lcmp          cmp for lambda derivative
c     lcphi         cphi for lambda derivative
c     lfmp          fmp for lambda derivative
c     lfphi         fphi for lambda derivative
c     poleorig      original pole
c     lqgrid        qgrid for lambda derivative
c     lmdatrial     flag to evaluate lambda bias for a trial move
c     use_abf       flag to use adaptive biasing force
c     use_abfdyn    flag to propagate abf lambda particle
c     use_dlmda     logical flag governing use of lambda derivative
c     use_edlmda    flag that the multipole term has a lambda deriv
c     use_elmdamap  flag that elambda follows the main lambda map
c     use_emdt      flag governing use of multipole dual topology
c     use_epdt      flag governing use of polarization dual topology
c     use_evdt      flag governing use of van der Waals dual topology
c     use_mainlmda  flag that a main lambda value was specified
c     use_meta      flag to use metadynamics
c     use_metadyn   flag to propagate metadynamics lambda particle
c     use_ost       flag to use orthogonal space tempering
c     use_ostdyn    flag to propagate lambda particle
c     use_pdlmda    flag that the polarization term has a lambda deriv
c     use_plmda     flag governing rescale to a decoupled plambda
c     use_plmdamap  flag that plambda follows the main lambda map
c     use_relstage  flag to use staged relative free energy schedule
c     use_ti        flag to use thermodynamic integration
c     use_vdlmda    flag that the van der Waals term has a lambda deriv
c     use_vlmdamap  flag that vlambda follows the main lambda map
c     douindorig    original douind
c     elmdamap      mapping type from main to electrostatic lambda
c     plmdamap      mapping type from main to polarization lambda
c     vlmdamap      mapping type from main to van der Waals lambda
c     lmdaengymode  free energy being computed, ABS or REL
c     lmdasampmode  lambda sampling method, OST, META, TI, ABF or NONE
c     relstage      declared leg of the staged schedule
c     abflabel      history label record of the abf history file
c     abftitle      title record of the abf history file
c     lmdasavefile  name of the file holding the lambda bias history
c
c
      module dlmda
      implicit none
      integer elmdaapmn
      integer elmdaexp
      integer elmdainvn
      integer emdtexp
      integer epdtexp
      integer erelst0
      integer erelst1
      integer evdtexp
      integer lmdaintv
      integer lmdanpa
      integer lmdanpb
      integer lmdanpc
      integer lmdastep
      integer nlmda
      integer nlmdahist
      integer nlmdasave
      integer nrelsub
      integer plmdaapmn
      integer plmdaexp
      integer plmdainvn
      integer prelst0
      integer prelst1
      integer rellig1
      integer rellig2
      integer relnone
      integer sizelmdahist
      integer vlmdaapmn
      integer vlmdaexp
      integer vlmdainvn
      integer vrelst0
      integer vrelst1
      integer, allocatable :: lmdaihist(:)
      parameter (nrelsub=5)
      parameter (rellig1=1)
      parameter (rellig2=2)
      parameter (relnone=3)
      real*8 d2edl2
      real*8 d2eldlmda2
      real*8 d2emdl2
      real*8 d2epdl2
      real*8 d2evdl2
      real*8 d2pldlmda2
      real*8 d2vldlmda2
      real*8 dedl
      real*8 dedlavg
      real*8 dedlstd
      real*8 deffdl
      real*8 deldlmda
      real*8 demdl
      real*8 depdl
      real*8 devdl
      real*8 dpldlmda
      real*8 dvldlmda
      real*8 elmdaapmrho
      real*8 elmdainveps
      real*8 lmdaavg
      real*8 lmdaddgdl
      real*8 lmdadeltag
      real*8 lmdadfdl
      real*8 lmdadt
      real*8 lmdafric
      real*8 lmdamass
      real*8 lmdaparatio
      real*8 lmdapbratio
      real*8 lmdapcratio
      real*8 lmdastd
      real*8 lmdatheta
      real*8 lmdavbias
      real*8 lmdavtheta
      real*8 plmdaapmrho
      real*8 plmdainveps
      real*8 qntelmda0
      real*8 qntelmda1
      real*8 qntplmda0
      real*8 qntplmda1
      real*8 qntvlmda0
      real*8 qntvlmda1
      real*8 vlmdaapmrho
      real*8 vlmdainveps
      real*8 wlmda
      real*8 wlmda2
      real*8 demvirdl(3,3)
      real*8 depvirdl(3,3)
      real*8 devvirdl(3,3)
      real*8 dvirdl(3,3)
      real*8, allocatable :: bdplorig(:)
      real*8, allocatable :: bflxorig(:)
      real*8, allocatable :: lmdafhist(:)
      real*8, allocatable :: lmdaflist(:)
      real*8, allocatable :: lmdafmean(:)
      real*8, allocatable :: lmdafsum(:)
      real*8, allocatable :: lmdafwt(:)
      real*8, allocatable :: lmdalhist(:)
      real*8, allocatable :: lmdallist(:)
      real*8, allocatable :: pchg0orig(:)
      real*8, allocatable :: pchgorig(:)
      real*8, allocatable :: pcoreorig(:)
      real*8, allocatable :: polarityorig(:)
      real*8, allocatable :: pval0orig(:)
      real*8, allocatable :: pvalorig(:)
      real*8, allocatable :: abflxorig(:,:)
      real*8, allocatable :: aflxorig(:,:)
      real*8, allocatable :: dfmdl(:,:)
      real*8, allocatable :: dfpdl(:,:)
      real*8, allocatable :: dfsumdl(:,:)
      real*8, allocatable :: dfvdl(:,:)
      real*8, allocatable :: lcmp(:,:)
      real*8, allocatable :: lcphi(:,:)
      real*8, allocatable :: lfmp(:,:)
      real*8, allocatable :: lfphi(:,:)
      real*8, allocatable :: poleorig(:,:)
      real*8, allocatable :: lqgrid(:,:,:,:)
      logical lmdatrial
      logical use_abf
      logical use_abfdyn
      logical use_dlmda
      logical use_edlmda
      logical use_elmdamap
      logical use_emdt
      logical use_epdt
      logical use_evdt
      logical use_mainlmda
      logical use_meta
      logical use_metadyn
      logical use_ost
      logical use_ostdyn
      logical use_pdlmda
      logical use_plmda
      logical use_plmdamap
      logical use_relstage
      logical use_ti
      logical use_vdlmda
      logical use_vlmdamap
      logical, allocatable :: douindorig(:)
      character*3 elmdamap
      character*3 plmdamap
      character*3 vlmdamap
      character*4 lmdaengymode
      character*4 lmdasampmode
      character*4 relstage
      character*40 abflabel
      character*40 abftitle
      character*240 lmdasavefile
      save
      end
