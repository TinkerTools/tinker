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
c     fli0           index of flambda bin where flambda = 0
c     nflmda         number of flambda bins
c     nmetahist      total number of metadynamics gaussians
c     nmethistsave   number of metadynamics gaussians written to file
c     sizemetahist   current allocation size for metadynamics gaussians
c     metaihist      lmdastep step at which each gaussian was added
c     osthist        packed lambda/flambda bin for saved gaussians
c     ostnext        next histogram index in the same lambda/flambda bin
c     osthead        first histogram index for each lambda/flambda bin
c     hbias          height of biasing gaussian
c     maxwfhist      maximum flambda width of histogram gaussians
c     maxwlhist      maximum lambda width of histogram gaussians
c     ostbdgdfl      saved bias dg/dflambda from eostbias for eostdyn
c     ostbdgdl       saved bias dg/dlambda from eostbias for eostdyn
c     ostcvdif       max drift between first and last convergence bin
c     ostcvrat       max ratio of sample deviation to sample average
c     ostcvslp       max fitted sample slope over a deposit interval
c     ostcvstd       max sample deviation over a deposit interval
c     ostdedlslp     fitted dU/dL change per sample over an interval
c     ostdgdl        current dg/dlambda value
c     ostgtempgamma  global tempering factor scaling height decay by kT
c     ostgthresh     global bias threshold for untempered heights
c     ostlambdaslp   fitted lambda change per sample over an interval
c     ostltempgamma  local tempering factor scaling height decay by kT
c     ostlthresh     local bias threshold for untempered heights
c     oststdev       gaussian cutoff distance in standard deviations
c     wfhist         flambda width of new histogram gaussians
c     wflmda         width of flambda bins
c     wflmda2        half width of flambda bins
c     wlhist         lambda width of new histogram gaussians
c     dvmetagrid     d(vmetagrid)/dlambda at each lambda bin
c     metahhist      height of metadynamics gaussians
c     metalhist      lambda center of metadynamics gaussians
c     metawhist      lambda width of metadynamics gaussians
c     osthhist       height of gaussians saved in histogram
c     ostwfhist      flambda width of gaussians saved in histogram
c     ostwlhist      lambda width of gaussians saved in histogram
c     vkernelmax     maximum gkernel over flambda at each lambda bin
c     vmetagrid      metadynamics bias at each lambda bin center
c     gfkernel       d(gkernel)/dflambda values on grid
c     gkernel        ost bias potential on the lambda/flambda grid
c     glfkernel      mixed derivative of gkernel on grid
c     glkernel       d(gkernel)/dlambda values on grid
c     fastkernel     flag to use fused g and f kernel updates
c     ostinterpol    flag to interpolate ost g kernel from grid
c     use_ostgtemp   flag to temper heights by the global bias level
c     use_ostltemp   flag to temper heights by the lambda bin excess
c     metasavefile   name of the file holding the metadynamics history
c     ostlabel       history label record of the ost history file
c     osttitle       title record of the ost history file
c
c
      module ost
      implicit none
      integer fli0
      integer nflmda
      integer nmetahist
      integer nmethistsave
      integer sizemetahist
      integer, allocatable :: metaihist(:)
      integer, allocatable :: osthist(:)
      integer, allocatable :: ostnext(:)
      integer, allocatable :: osthead(:,:)
      real*8 hbias
      real*8 maxwfhist
      real*8 maxwlhist
      real*8 ostbdgdfl
      real*8 ostbdgdl
      real*8 ostcvdif
      real*8 ostcvrat
      real*8 ostcvslp
      real*8 ostcvstd
      real*8 ostdedlslp
      real*8 ostdgdl
      real*8 ostgtempgamma
      real*8 ostgthresh
      real*8 ostlambdaslp
      real*8 ostltempgamma
      real*8 ostlthresh
      real*8 oststdev
      real*8 wfhist
      real*8 wflmda
      real*8 wflmda2
      real*8 wlhist
      real*8, allocatable :: dvmetagrid(:)
      real*8, allocatable :: metahhist(:)
      real*8, allocatable :: metalhist(:)
      real*8, allocatable :: metawhist(:)
      real*8, allocatable :: osthhist(:)
      real*8, allocatable :: ostwfhist(:)
      real*8, allocatable :: ostwlhist(:)
      real*8, allocatable :: vkernelmax(:)
      real*8, allocatable :: vmetagrid(:)
      real*8, allocatable :: gfkernel(:,:)
      real*8, allocatable :: gkernel(:,:)
      real*8, allocatable :: glfkernel(:,:)
      real*8, allocatable :: glkernel(:,:)
      logical fastkernel
      logical ostinterpol
      logical use_ostgtemp
      logical use_ostltemp
      character*40 ostlabel
      character*40 osttitle
      character*240 metasavefile
      parameter (ostlabel=' Gaussian History :')
      parameter (osttitle=' Orthogonal Space Tempering History :')
      save
      end
