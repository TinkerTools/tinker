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
c     nmethistsave  number of metadynamics gaussians written to restart
c     nosthist      total number of histograms
c     nosthistsave  number of histograms written to the restart file
c     osteqratio    fraction of hist interval to equilibrate
c     ostcvbin      convergence sub-bins per gaussian deposit interval
c     ostnavg       samples averaged between hist updates
c     ostnequil     samples skipped before hist averaging
c     sizeosthist   current allocation size for saved histograms
c     sizemetahist  current allocation size for metadynamics gaussians
c     metaihist     iost step at which each gaussian was added
c     osthead       first histogram index for each lambda/flambda bin
c     osthist       packed lambda/flambda bin for saved gaussians
c     ostihist      iost step at which each gaussian was added
c     ostnext       next histogram index in the same lambda/flambda bin
c     deffdl        effective lambda derivative for propagation
c     eosttot       total ost free energy
c     hbias         height of biasing gaussian
c     maxwfhist     maximum flambda width of histogram gaussians
c     maxwlhist     maximum lambda width of histogram gaussians
c     ostbdfdl      saved bias df/dlambda from eostbias for eostdyn
c     ostbdgdfl     saved bias dg/dflambda from eostbias for eostdyn
c     ostbdgdl      saved bias dg/dlambda from eostbias for eostdyn
c     ostbvbias     saved bias energy shift from eostbias
c     ostcvdif      max drift between first and last convergence bin
c     ostcvrat      max ratio of sample deviation to sample average
c     ostcvslp      max fitted sample slope over a deposit interval
c     ostcvstd      max sample deviation over a deposit interval
c     ostddgdl      current dDeltaG/dlambda value
c     ostdedl       current unbiased dU/dlambda value
c     ostdedlavg    average dE/dlambda value between hist updates
c     ostdedlstd    standard deviation of dE/dlambda between updates
c     ostdgdl       current dg/dlambda value
c     ostdt         time step for theta lambda propagation
c     ostfriction   friction coefficient for theta lambda coordinate
c     ostlambda     main lambda value in orthogonal space sampling
c     ostlambdaavg  average main lambda value between hist updates
c     ostlambdastd  standard deviation of lambda between hist updates
c     ostmass       fictitious mass of theta lambda coordinate
c     oststdev      gaussian cutoff distance in standard deviations
c     osttheta      theta coordinate used to propagate lambda
c     ostvtheta     velocity of the theta lambda coordinate
c     tempergamma   tempering factor scaling the height decay by kT
c     temperthresh  bias level below which heights are untempered
c     wfhist        flambda width of new histogram gaussians
c     wflmda        width of flambda bins
c     wflmda2       half width of flambda bins
c     wlhist        lambda width of new histogram gaussians
c     wlmda         width of lambda bins
c     wlmda2        half width of lambda bins
c     fkernel       free energy mean force at each lambda bin
c     fsumkernel    numerator of free energy mean force kernel
c     gfkernel      d(gkernel)/dflambda values on grid
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
c     glfkernel     mixed derivative of gkernel on grid
c     glkernel      d(gkernel)/dlambda values on grid
c     pfkernel      partition function for free energy mean force
c     fastkernel    flag to use fused g and f kernel updates
c     metarestart   flag to indicate metadynamics restart data was read
c     ostemper      flag to temper the deposited gaussian heights
c     ostinterpol   flag to interpolate ost g kernel from grid
c     ostrestart    flag to indicate ost restart data was read
c     osttrial      flag to evaluate ost bias w/o depositing gaussians
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
      integer nmethistsave
      integer nosthist
      integer nosthistsave
      integer ostcvbin
      integer ostnavg
      integer ostnequil
      integer sizeosthist
      integer sizemetahist
      integer, allocatable :: metaihist(:)
      integer, allocatable :: osthead(:,:)
      integer, allocatable :: osthist(:)
      integer, allocatable :: ostihist(:)
      integer, allocatable :: ostnext(:)
      real*8 deffdl
      real*8 eosttot
      real*8 hbias
      real*8 maxwfhist
      real*8 maxwlhist
      real*8 ostbdfdl
      real*8 ostbdgdfl
      real*8 ostbdgdl
      real*8 ostbvbias
      real*8 ostcvdif
      real*8 ostcvrat
      real*8 ostcvslp
      real*8 ostcvstd
      real*8 ostddgdl
      real*8 ostdedl
      real*8 ostdedlavg
      real*8 ostdedlstd
      real*8 ostdgdl
      real*8 ostdt
      real*8 osteqratio
      real*8 ostfriction
      real*8 ostlambda
      real*8 ostlambdaavg
      real*8 ostlambdastd
      real*8 ostmass
      real*8 oststdev
      real*8 osttheta
      real*8 ostvtheta
      real*8 tempergamma
      real*8 temperthresh
      real*8 wfhist
      real*8 wflmda
      real*8 wflmda2
      real*8 wlhist
      real*8 wlmda
      real*8 wlmda2
      real*8, allocatable :: fkernel(:)
      real*8, allocatable :: fsumkernel(:)
      real*8, allocatable :: gfkernel(:,:)
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
      real*8, allocatable :: glfkernel(:,:)
      real*8, allocatable :: glkernel(:,:)
      real*8, allocatable :: pfkernel(:)
      logical fastkernel
      logical metarestart
      logical ostemper
      logical ostinterpol
      logical ostrestart
      logical osttrial
      save
      end
