c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  module ost  --  orthogonal space tempering variables  ##
c     ##                                                        ##
c     ############################################################
c
c
c     fli0          index of flambda bin where flambda = 0
c     iost          current md step while running ost (=istep)
c     iosthist      steps between histogram additions
c     nflmda        number of flambda bins
c     nlmda         number of lambda bins
c     nmetahist     total number of metadynamics gaussians
c     nosthist      total number of histograms
c     osteqratio    fraction of hist interval to equilibrate
c     ostnavg       samples averaged between hist updates
c     ostnequil     samples skipped before hist averaging
c     ostemexp      exponent for electrostatic exponential mapping
c     ostepexp      exponent for polarization exponential mapping
c     ostevexp      exponent for van der Waals exponential mapping
c     ostinvemn     inverse-power exponent for electrostatic mapping
c     ostinvepn     inverse-power exponent for polarization mapping
c     ostinvevn     inverse-power exponent for van der Waals mapping
c     sizeosthist   current allocation size for saved histograms
c     sizemetahist  current allocation size for metadynamics gaussians
c     osthist       packed lambda/flambda bin for saved gaussians
c     ostnext       next histogram index in the same lambda/flambda bin
c     osthead       first histogram index for each lambda/flambda bin
c     deffdl        effective lambda derivative for propagation
c     eosttot       total ost free energy
c     hbias         height of biasing gaussian
c     maxwfhist     maximum flambda width of histogram gaussians
c     maxwlhist     maximum lambda width of histogram gaussians
c     ostddgdl      current dDeltaG/dlambda value
c     ostdedl       current unbiased dU/dlambda value
c     ostdedlavg    average dE/dlambda value between hist updates
c     ostdedlstd    standard deviation of dE/dlambda between updates
c     ostdgdl       current dg/dlambda value
c     ostdt         time step for theta lambda propagation
c     ostelmda0     sublambda lower bound for electrostatics
c     ostelmda1     sublambda upper bound for electrostatics
c     ostfriction   friction coefficient for theta lambda coordinate
c     ostinvemeps   shift for electrostatic inverse-power mapping
c     ostinvepeps   shift for polarization inverse-power mapping
c     ostinveveps   shift for van der Waals inverse-power mapping
c     ostlambda     main lambda value in orthogonal space sampling
c     ostlambdaavg  average main lambda value between hist updates
c     ostlambdastd  standard deviation of lambda between hist updates
c     ostmass       fictitious mass of theta lambda coordinate
c     ostplmda0     sublambda lower bound for polarization
c     ostplmda1     sublambda upper bound for polarization
c     oststdev      gaussian cutoff distance in standard deviations
c     osttheta      theta coordinate used to propagate lambda
c     ostvlmda0     sublambda lower bound for van der Waals
c     ostvlmda1     sublambda upper bound for van der Waals
c     ostvtheta     velocity of the theta lambda coordinate
c     wfhist        flambda width of new histogram gaussians
c     wflmda        width of flambda bins
c     wflmda2       half width of flambda bins
c     wlhist        lambda width of new histogram gaussians
c     wlmda         width of lambda bins
c     wlmda2        half width of lambda bins
c     fkernel       free energy mean force at each lambda bin
c     metahhist     height of metadynamics gaussians
c     metalhist     lambda center of metadynamics gaussians
c     metawhist     lambda width of metadynamics gaussians
c     ostfhist      flambda center of gaussians saved in histogram
c     osthhist      height of gaussians saved in histogram
c     ostflist      dE/dlambda values saved between hist updates
c     ostlhist      lambda center of gaussians saved in histogram
c     ostllist      lambda values saved between hist updates
c     ostwfhist     flambda width of gaussians saved in histogram
c     ostwlhist     lambda width of gaussians saved in histogram
c     gkernel       ost bias potential on the lambda/flambda grid
c     metarestart   flag to indicate metadynamics restart data was read
c     ostrestart    flag to indicate ost restart data was read
c     use_meta      flag to use metadynamics
c     use_metadyn   flag to propagate metadynamics lambda particle
c     use_ost       flag to use orthogonal space tempering
c     use_ostdyn    flag to propagate lambda particle
c     use_pol4f     flag to compute polarization lambda deriv for lmda=1
c     use_pol4i     flag to compute polarization lambda deriv for lmda=0
c     ostemap       mapping type from ost lambda to electrostatic lambda
c     ostpmap       mapping type from ost lambda to polarization lambda
c     ostvmap       mapping type from ost lambda to van der Waals lambda
c
c
      module ost
      implicit none
      integer fli0
      integer iost
      integer iosthist
      integer nflmda
      integer nlmda
      integer nmetahist
      integer nosthist
      integer ostnavg
      integer ostnequil
      integer ostemexp
      integer ostepexp
      integer ostevexp
      integer ostinvemn
      integer ostinvepn
      integer ostinvevn
      integer sizeosthist
      integer sizemetahist
      integer, allocatable :: osthist(:)
      integer, allocatable :: ostnext(:)
      integer, allocatable :: osthead(:,:)
      real*8 deffdl
      real*8 eosttot
      real*8 hbias
      real*8 maxwfhist
      real*8 maxwlhist
      real*8 ostddgdl
      real*8 ostdedl
      real*8 ostdedlavg
      real*8 ostdedlstd
      real*8 ostdgdl
      real*8 ostdt
      real*8 ostelmda0
      real*8 ostelmda1
      real*8 osteqratio
      real*8 ostfriction
      real*8 ostinvemeps
      real*8 ostinvepeps
      real*8 ostinveveps
      real*8 ostlambda
      real*8 ostlambdaavg
      real*8 ostlambdastd
      real*8 ostmass
      real*8 ostplmda0
      real*8 ostplmda1
      real*8 oststdev
      real*8 osttheta
      real*8 ostvlmda0
      real*8 ostvlmda1
      real*8 ostvtheta
      real*8 wfhist
      real*8 wflmda
      real*8 wflmda2
      real*8 wlhist
      real*8 wlmda
      real*8 wlmda2
      real*8, allocatable :: fkernel(:)
      real*8, allocatable :: metahhist(:)
      real*8, allocatable :: metalhist(:)
      real*8, allocatable :: metawhist(:)
      real*8, allocatable :: ostflist(:)
      real*8, allocatable :: ostfhist(:)
      real*8, allocatable :: osthhist(:)
      real*8, allocatable :: ostllist(:)
      real*8, allocatable :: ostlhist(:)
      real*8, allocatable :: ostwfhist(:)
      real*8, allocatable :: ostwlhist(:)
      real*8, allocatable :: gkernel(:,:)
      logical metarestart
      logical ostrestart
      logical use_meta
      logical use_metadyn
      logical use_ost
      logical use_ostdyn
      logical use_pol4f
      logical use_pol4i
      character*3 ostemap
      character*3 ostpmap
      character*3 ostvmap
      save
      end
