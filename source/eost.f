c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine eostbias -- apply ost/metadynamics bias  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "eostbias" evaluates the orthogonal space tempering or
c     metadynamics bias at the current configuration and folds it
c     into the energy, Cartesian gradient, and virial.  It has no
c     side effects on the histogram or the lambda particle, so it
c     may be called from a Monte Carlo barostat trial as well as
c     from a normal dynamics step.  The bias derivatives needed to
c     propagate the lambda particle are saved for a later elmdadyn
c     call in the same step.
c
c
      subroutine eostbias
      use atoms
      use deriv
      use dlmda
      use energi
      use mutant
      use ost
      use virial
      implicit none
      integer i,j
      real*8 egbias,dgdl,dgdfl
      real*8 eostlmda,dfdl
      real*8 vbias,dvdl
c
c
c     the metadynamics bias depends on lambda alone
c
      if (use_metadyn) then
         call emetabias (lambda,vbias,dvdl)
         esum = esum + vbias
         lmdavbias = vbias
         ostbdgdl = dvdl
         return
      end if
c
c     evaluate the ost g kernel bias and free energy f term
c
      if (ostinterpol) then
         call egkernelinterpolate (egbias,dgdl,dgdfl)
      else
         call egkernel (egbias,dgdl,dgdfl)
      end if
      call efreelmda (lambda,eostlmda,dfdl)
      esum = esum + egbias - eostlmda
      lmdavbias = egbias - eostlmda
      ostbdgdl = dgdl
      ostbdgdfl = dgdfl
      lmdadfdl = dfdl
c
c     add the second order lambda force and virial from the bias
c
      do i = 1, n
         do j = 1, 3
            desum(j,i) = desum(j,i) + dgdfl*dfsumdl(j,i)
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = vir(j,i) + dgdfl*dvirdl(j,i)
         end do
      end do
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine ostdeposit -- add one ost bias gaussian  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "ostdeposit" records one accepted interval of orthogonal space
c     tempering as a biasing gaussian centered on the interval lambda
c     and dU/dlambda averages, then updates the bias kernels and the
c     free energy estimate
c
c
      subroutine ostdeposit
      use dlmda
      use ost
      implicit none
      integer k
      integer ilmda,iflmda
      integer lmdabin,flambdabin
      real*8 efreetot
      real*8 ostvminimax
      real*8 temperedheight
c
c
c     ensure histogram contains the unbiased dU/dlambda value
c
      ilmda = lmdabin(lmdaavg)
      maxwlhist = max(maxwlhist,wlhist)
      maxwfhist = max(maxwfhist,wfhist)
      call ensureflambda (dedlavg)
      iflmda = flambdabin(dedlavg)
c
c     ensure histogram array is sufficiently large
c
      nlmdahist = nlmdahist + 1
      if (nlmdahist .gt. sizelmdahist)  call resizeosthist
      call ij_to_k(ilmda,iflmda,nlmda,k)
c
c     save histogram information
c
      osthist(nlmdahist) = k
      lmdaihist(nlmdahist) = lmdastep
      lmdalhist(nlmdahist) = lmdaavg
      lmdafhist(nlmdahist) = dedlavg
      osthhist(nlmdahist) = temperedheight(ostvminimax(),
     &                                    vkernelmax(ilmda))
      ostwlhist(nlmdahist) = wlhist
      ostwfhist(nlmdahist) = wfhist
      ostnext(nlmdahist) = osthead(ilmda,iflmda)
      osthead(ilmda,iflmda) = nlmdahist
      if (fastkernel) then
         call updatekernels
      else
         call updategkernel
         call buildfkernel
      end if
      lmdadeltag = efreetot()
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine metadeposit -- add one metadynamics gaussian  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "metadeposit" records one interval of lambda metadynamics as a
c     gaussian centered on the interval lambda average, then updates
c     the bin lookup, the bias grid and the free energy estimate
c
c     note only the global tempering factor is applied, since the
c     free energy is read from the bias itself, and a local factor
c     would keep the bias from filling barriers above its threshold
c
c
      subroutine metadeposit
      use dlmda
      use ost
      implicit none
      integer ilmda
      integer lmdabin
      real*8 vglobal
      real*8 metadeltag
      real*8 metavminimax
      real*8 temperedheight
c
c
c     save the new gaussian with the global tempering factor alone
c
      nmetahist = nmetahist + 1
      if (nmetahist .gt. sizemetahist)  call resizemeta
      vglobal = metavminimax ()
      metalhist(nmetahist) = lmdaavg
      metahhist(nmetahist) = temperedheight(vglobal,vglobal)
      metawhist(nmetahist) = wlhist
      metaihist(nmetahist) = lmdastep
      maxwlhist = max(maxwlhist,wlhist)
c
c     link the gaussian into its lambda bin, then add it to the grid
c
      ilmda = lmdabin(lmdaavg)
      metanext(nmetahist) = metahead(ilmda)
      metahead(ilmda) = nmetahist
      call addmetagrid (nmetahist)
      lmdadeltag = metadeltag()
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  function ostvminimax -- global ost path bias level  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "ostvminimax" returns the minimum over lambda of the maximum
c     ost bias reached along flambda, the bias level that any path
c     across the lambda range has to have been filled to
c
c
      function ostvminimax ()
      use dlmda
      use ost
      implicit none
      integer i
      real*8 ostvminimax
c
c
c     reduce the running per lambda bin maxima
c
      ostvminimax = 0.0d0
      if (nlmda .lt. 1)  return
      if (.not. allocated(vkernelmax))  return
      ostvminimax = vkernelmax(1)
      do i = 2, nlmda
         ostvminimax = min(ostvminimax,vkernelmax(i))
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function metavminimax -- global metadynamics bias level  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "metavminimax" returns the minimum over lambda of the deposited
c     metadynamics bias, the level any path across lambda is filled to
c
c
      function metavminimax ()
      use dlmda
      use ost
      implicit none
      integer i
      real*8 metavminimax
c
c
c     reduce the metadynamics bias over the lambda bin centers
c
      metavminimax = 0.0d0
      if (nlmda .lt. 1)  return
      if (.not. allocated(vmetagrid))  return
      metavminimax = vmetagrid(1)
      do i = 2, nlmda
         metavminimax = min(metavminimax,vmetagrid(i))
      end do
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  function temperedheight -- tempered gaussian height  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "temperedheight" scales down the height of a new biasing
c     gaussian by a global factor once the path bias level passes
c     its threshold, and by a local factor once the bias level at
c     the deposit lambda bin exceeds the path bias level by more
c     than the local threshold, each decaying on a scale of kT
c     times its own tempering factor
c
c
      function temperedheight (vglobal,vlocal)
      use bath
      use ost
      use units
      implicit none
      real*8 temperedheight
      real*8 vglobal,vlocal
      real*8 rt,expo,excess
c
c
c     an untempered run always deposits at the full height
c
      temperedheight = hbias
      rt = gasconst * kelvin
      if (rt .le. 0.0d0)  return
      expo = 0.0d0
c
c     the global factor follows the path bias level
c
      if (use_ostgtemp .and. ostgtempgamma.gt.0.0d0) then
         excess = max(0.0d0,vglobal-ostgthresh)
         expo = expo - excess/(rt*ostgtempgamma)
      end if
c
c     the local factor follows the excess of the deposit bin
c
      if (use_ostltemp .and. ostltempgamma.gt.0.0d0) then
         excess = max(0.0d0,vlocal-vglobal-ostlthresh)
         expo = expo - excess/(rt*ostltempgamma)
      end if
      temperedheight = hbias * exp(expo)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine emetabias -- evaluate metadynamics bias  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "emetabias" evaluates Vbias(lambda) and dVbias/dlambda for the
c     sum of one-dimensional normalized metadynamics gaussians, each
c     counted along with its two images reflected at the lambda end
c     points so that the bias has no gradient through the walls
c
c     note only gaussians in the lambda bins within the cutoff of
c     "oststdev" widths are visited, through the bin lookup, and only
c     the first "nmetahist" gaussians are counted, so a history being
c     replayed gaussian by gaussian gives the bias at each stage
c
c
      subroutine emetabias (lambda,vbias,dvdl)
      use dlmda
      use math
      use ost
      implicit none
      integer ihist
      integer img,nimg
      integer ilmda,klmda
      integer lcenter,lcount
      integer nlcut
      integer lmdabin
      integer nbincut
      real*8 lambda
      real*8 vbias,dvdl
      real*8 delta
      real*8 sig,sig2
      real*8 bias,pref
      real*8 source
      real*8 lmdaimgsrc
c
c
c     initialize metadynamics bias and derivative
c
      vbias = 0.0d0
      dvdl = 0.0d0
c
c     use the accumulated grid when interpolation was requested
c
      if (ostinterpol .and. nmetahist.gt.0) then
         call emetabiasinterpolate (lambda,vbias,dvdl)
         return
      end if
c
c     use the largest gaussian width for a conservative bin cutoff
c
      ilmda = lmdabin(lambda)
      nlcut = nbincut(oststdev*maxwlhist,wlmda) + 1
c
c     loop over nearby lambda bins and their mirror images
c
      do klmda = -nlcut, nlcut
         lcenter = ilmda + klmda
         call lmdaimgbin (lcenter,lcount,nimg)
c
c     add each gaussian of the bin with the image it stands for
c
         if (lcount.ge.1 .and. lcount.le.nlmda) then
            ihist = metahead(lcount)
            do while (ihist .ne. 0)
               sig = metawhist(ihist)
               if (ihist.le.nmetahist .and. sig.gt.0.0d0) then
                  sig2 = sig * sig
                  pref = metahhist(ihist) / (sig*sqrt(2.0d0*pi))
                  do img = 1, nimg
                     source = lmdaimgsrc(lcenter,img,metalhist(ihist))
                     delta = lambda - source
                     if (abs(delta) .le. oststdev*sig) then
                        bias = pref * exp(-0.5d0*delta*delta/sig2)
                        vbias = vbias + bias
                        dvdl = dvdl - delta*bias/sig2
                     end if
                  end do
               end if
               ihist = metanext(ihist)
            end do
         end if
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emetabiasinterpolate -- interpolated meta bias  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emetabiasinterpolate" evaluates the metadynamics bias and its
c     lambda derivative by cubic Hermite interpolation of the values
c     accumulated at the lambda bin centers
c
c
      subroutine emetabiasinterpolate (lambda,vbias,dvdl)
      use dlmda
      use ost
      implicit none
      integer i,ia
      integer il0
      real*8 lambda
      real*8 vbias,dvdl
      real*8 lam,l0,x
      real*8 x2,x3
      real*8 val,der
      real*8 hxv(2),hxd(2)
      real*8 dhxv(2),dhxd(2)
c
c
c     initialize metadynamics bias and derivative
c
      vbias = 0.0d0
      dvdl = 0.0d0
      if (nlmda .lt. 2)  return
c
c     find the lambda interval holding the requested value
c
      lam = min(1.0d0,max(0.0d0,lambda))
      if (lam .ge. 1.0d0) then
         il0 = nlmda - 1
      else
         il0 = int(lam/wlmda) + 1
         il0 = max(1,min(il0,nlmda-1))
      end if
      l0 = dble(il0-1) * wlmda
      x = (lam-l0) / wlmda
c
c     cubic Hermite basis functions and their derivatives
c
      x2 = x * x
      x3 = x2 * x
      hxv(1) = 2.0d0*x3 - 3.0d0*x2 + 1.0d0
      hxv(2) = -2.0d0*x3 + 3.0d0*x2
      hxd(1) = x3 - 2.0d0*x2 + x
      hxd(2) = x3 - x2
      dhxv(1) = 6.0d0*x2 - 6.0d0*x
      dhxv(2) = -6.0d0*x2 + 6.0d0*x
      dhxd(1) = 3.0d0*x2 - 4.0d0*x + 1.0d0
      dhxd(2) = 3.0d0*x2 - 2.0d0*x
c
c     accumulate the value and slope contributions of both nodes
c
      do ia = 1, 2
         i = il0 + ia - 1
         val = vmetagrid(i)
         der = wlmda * dvmetagrid(i)
         vbias = vbias + hxv(ia)*val + hxd(ia)*der
         dvdl = dvdl + dhxv(ia)*val + dhxd(ia)*der
      end do
      dvdl = dvdl / wlmda
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine addmetagrid -- accumulate one meta gaussian  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "addmetagrid" adds one deposited metadynamics gaussian and its
c     two reflected images to the bias and derivative values stored
c     at the lambda bin centers within the cutoff of "oststdev" widths
c
c
      subroutine addmetagrid (ihist)
      use dlmda
      use math
      use ost
      implicit none
      integer ihist
      integer il,img
      integer il1,il2
      real*8 lambda,delta
      real*8 sig,sig2
      real*8 bias,pref
      real*8 cut,x1,x2
      real*8 src(3)
c
c
c     a zero width gaussian makes no contribution
c
      sig = metawhist(ihist)
      if (sig .le. 0.0d0)  return
      sig2 = sig * sig
      pref = metahhist(ihist) / (sig*sqrt(2.0d0*pi))
      cut = oststdev * sig
      src(1) = metalhist(ihist)
      src(2) = -metalhist(ihist)
      src(3) = 2.0d0 - metalhist(ihist)
c
c     spread the gaussian and its images over the bin centers
c     within the cutoff
c
      do img = 1, 3
         x1 = (src(img)-cut)/wlmda + 1.0d0
         x2 = (src(img)+cut)/wlmda + 1.0d0
         if (x1.le.dble(nlmda) .and. x2.ge.1.0d0) then
            il1 = int(max(1.0d0,x1))
            il2 = min(nlmda,int(min(dble(nlmda),x2))+1)
            do il = il1, il2
               lambda = dble(il-1) * wlmda
               delta = lambda - src(img)
               if (abs(delta) .le. cut) then
                  bias = pref * exp(-0.5d0*delta*delta/sig2)
                  vmetagrid(il) = vmetagrid(il) + bias
                  dvmetagrid(il) = dvmetagrid(il) - delta*bias/sig2
               end if
            end do
         end if
      end do
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine buildmetaindex -- build metadynamics bin lookup  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "buildmetaindex" rebuilds the linked-list lookup of the saved
c     metadynamics gaussians by lambda bin, in the order deposited
c
c
      subroutine buildmetaindex
      use dlmda
      use ost
      implicit none
      integer i
      integer ihist
      integer ilmda
      integer lmdabin
c
c
c     clear the bin heads and next links
c
      do i = 1, nlmda
         metahead(i) = 0
      end do
      do i = 1, sizemetahist
         metanext(i) = 0
      end do
c
c     insert each saved gaussian into its lambda bin
c
      do ihist = 1, nmetahist
         ilmda = lmdabin(metalhist(ihist))
         metanext(ihist) = metahead(ilmda)
         metahead(ilmda) = ihist
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function metadeltag -- metadynamics free energy estimate  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "metadeltag" returns G(1)-G(0) = -Vbias(1) + Vbias(0)
c
c     note the end point bias is read from the grid, which holds the
c     exact bias at the bin centers for every gaussian added so far
c
c
      function metadeltag ()
      use dlmda
      use ost
      implicit none
      real*8 metadeltag
c
c
c     take the end point metadynamics bias values from the grid
c
      metadeltag = 0.0d0
      if (.not. allocated(vmetagrid))  return
      metadeltag = -vmetagrid(nlmda) + vmetagrid(1)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine resizemeta -- resize metadynamics history  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "resizemeta" doubles storage for saved metadynamics gaussians
c
c
      subroutine resizemeta
      use ost
      implicit none
      integer oldsize,newsize
      integer, allocatable :: itmp(:)
      real*8, allocatable :: rtmp(:)
c
c
c     grow each array, keeping the saved entries and zeroing the rest
c
      oldsize = sizemetahist
      newsize = 2 * oldsize
      allocate (rtmp(newsize))
      rtmp(1:oldsize) = metalhist(1:oldsize)
      rtmp(oldsize+1:newsize) = 0.0d0
      call move_alloc (rtmp,metalhist)
      allocate (rtmp(newsize))
      rtmp(1:oldsize) = metahhist(1:oldsize)
      rtmp(oldsize+1:newsize) = 0.0d0
      call move_alloc (rtmp,metahhist)
      allocate (rtmp(newsize))
      rtmp(1:oldsize) = metawhist(1:oldsize)
      rtmp(oldsize+1:newsize) = 0.0d0
      call move_alloc (rtmp,metawhist)
      allocate (itmp(newsize))
      itmp(1:oldsize) = metaihist(1:oldsize)
      itmp(oldsize+1:newsize) = 0
      call move_alloc (itmp,metaihist)
      allocate (itmp(newsize))
      itmp(1:oldsize) = metanext(1:oldsize)
      itmp(oldsize+1:newsize) = 0
      call move_alloc (itmp,metanext)
      sizemetahist = newsize
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine ensureflambda -- resize flambda histogram  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "ensureflambda" grows the flambda grid if the current
c     dE/dlambda value falls outside the saved histogram range
c
c
      subroutine ensureflambda (dudl)
      use dlmda
      use ost
      implicit none
      integer i
      integer iflmda
      integer naddlow
      integer naddhigh
      integer nneed
      integer nchunk
      integer nfcut
      integer nbuffer
      integer oldnflmda
      integer oldfli0
      integer offset
      integer jflmda
      real*8 dudl
      real*8, allocatable :: gfkernel0(:,:)
      real*8, allocatable :: gkernel0(:,:)
      real*8, allocatable :: glfkernel0(:,:)
      real*8, allocatable :: glkernel0(:,:)
c
c
c     get the raw flambda bin before any endpoint clamping
c
      iflmda = nint(dudl / wflmda) + fli0
c
c     include the full gaussian tail plus a buffer on each side
c
      nchunk = 100
      nbuffer = nchunk
      nfcut = int(oststdev*maxwfhist/wflmda) + 1 + nbuffer
      if (iflmda-nfcut .ge. 1 .and.
     &    iflmda+nfcut .le. nflmda)  return
c
c     increase the flambda grid in chunks of bins
c
      oldnflmda = nflmda
      oldfli0 = fli0
      naddlow = 0
      naddhigh = 0
      if (iflmda-nfcut .lt. 1) then
         nneed = 1 - (iflmda - nfcut)
         naddlow = ((nneed + nchunk - 1) / nchunk) * nchunk
         nflmda = nflmda + naddlow
         fli0 = fli0 + naddlow
         iflmda = iflmda + naddlow
      end if
      if (iflmda+nfcut .gt. nflmda) then
         nneed = iflmda + nfcut - nflmda
         naddhigh = ((nneed + nchunk - 1) / nchunk) * nchunk
         nflmda = nflmda + naddhigh
      end if
      offset = fli0 - oldfli0
c
c     resize and rebuild the lookup index for saved histograms
c
      deallocate (osthead)
      allocate (osthead(nlmda,nflmda))
      call buildostindex
c
c     resize the flambda-dependent kernels while preserving old data
c
      allocate (gfkernel0(nlmda,oldnflmda))
      allocate (gkernel0(nlmda,oldnflmda))
      allocate (glfkernel0(nlmda,oldnflmda))
      allocate (glkernel0(nlmda,oldnflmda))
      do i = 1, nlmda
         do jflmda = 1, oldnflmda
            gfkernel0(i,jflmda) = gfkernel(i,jflmda)
            gkernel0(i,jflmda) = gkernel(i,jflmda)
            glfkernel0(i,jflmda) = glfkernel(i,jflmda)
            glkernel0(i,jflmda) = glkernel(i,jflmda)
         end do
      end do
      deallocate (gfkernel)
      deallocate (gkernel)
      deallocate (glfkernel)
      deallocate (glkernel)
      allocate (gfkernel(nlmda,nflmda))
      allocate (gkernel(nlmda,nflmda))
      allocate (glfkernel(nlmda,nflmda))
      allocate (glkernel(nlmda,nflmda))
      do i = 1, nlmda
         do jflmda = 1, nflmda
            gfkernel(i,jflmda) = 0.0d0
            gkernel(i,jflmda) = 0.0d0
            glfkernel(i,jflmda) = 0.0d0
            glkernel(i,jflmda) = 0.0d0
         end do
      end do
      do i = 1, nlmda
         do jflmda = 1, oldnflmda
            gfkernel(i,jflmda+offset) = gfkernel0(i,jflmda)
            gkernel(i,jflmda+offset) = gkernel0(i,jflmda)
            glfkernel(i,jflmda+offset) = glfkernel0(i,jflmda)
            glkernel(i,jflmda+offset) = glkernel0(i,jflmda)
         end do
      end do
      deallocate (gfkernel0)
      deallocate (gkernel0)
      deallocate (glfkernel0)
      deallocate (glkernel0)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine egkernel -- ost g energy and derivatives  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "egkernel" computes the ost biasing g kernel energy and
c     derivatives at the current lambda and unbiased dU/dlambda values
c
c
      subroutine egkernel (egbias,dgdl,dgdfl)
      use dlmda
      use math
      use mutant
      use ost
      implicit none
      integer ihist
      integer ilmda,iflmda
      integer klmda,kflmda
      integer lcenter,lcount
      integer flcenter
      integer nlcut,nfcut
      integer img,nimg
      integer lmdabin
      integer nbincut
      real*8 egbias,dgdl,dgdfl
      real*8 sigl,sigf
      real*8 sigl2,sigf2
      real*8 sigl2inv,sigf2inv
      real*8 pref
      real*8 sourcel,sourcefl
      real*8 ldelta,ldelta2
      real*8 fldelta,fldelta2
      real*8 bias,expl,expfl
      real*8 lmdaimgsrc
      logical getexpfl
c
c
c     initialize bias and derivatives
c
      egbias = 0.0d0
      dgdl = 0.0d0
      dgdfl = 0.0d0
c
c     get the current lambda and flambda bin indices
c
      ilmda = lmdabin(lambda)
      iflmda = nint(dedl / wflmda) + fli0
      if (iflmda .lt. 1 .or. iflmda .gt. nflmda)  return
c
c     use max gaussian widths for conservative bin cutoffs
c
      nlcut = nbincut(oststdev*maxwlhist,wlmda) + 1
      nfcut = nbincut(oststdev*maxwfhist,wflmda) + 1
c
c     loop over nearby lambda bins and their mirror images
c
      do klmda = -nlcut, nlcut
         lcenter = ilmda + klmda
         call lmdaimgbin (lcenter,lcount,nimg)
c
c     loop over nearby flambda bins without mirror images
c
         if (lcount .ge. 1 .and. lcount .le. nlmda) then
            do kflmda = -nfcut, nfcut
               flcenter = iflmda + kflmda
               if (flcenter .ge. 1 .and. flcenter .le. nflmda) then
                  ihist = osthead(lcount,flcenter)
                  do while (ihist .ne. 0)
                     sigl = ostwlhist(ihist)
                     sigf = ostwfhist(ihist)
                     sourcefl = lmdafhist(ihist)
                     fldelta = dedl - sourcefl
                     if (abs(fldelta) .le. oststdev*sigf) then
                        sigl2 = sigl * sigl
                        sigf2 = sigf * sigf
                        sigl2inv = 1.0d0 / sigl2
                        sigf2inv = 1.0d0 / sigf2
                        pref = osthhist(ihist)
     &                            / (2.0d0*pi*sigl*sigf)
                        getexpfl = .true.
                        do img = 1, nimg
                           sourcel = lmdaimgsrc(lcenter,img,
     &                                          lmdalhist(ihist))
                           ldelta = lambda - sourcel
                           if (abs(ldelta) .le. oststdev*sigl) then
c
c     the flambda factor is found once, when an image first counts
c
                              if (getexpfl) then
                                 fldelta2 = fldelta * fldelta
                                 expfl = exp(-0.5d0*fldelta2*sigf2inv)
                                 getexpfl = .false.
                              end if
                              ldelta2 = ldelta * ldelta
                              expl = exp(-0.5d0*ldelta2*sigl2inv)
                              bias = pref * expl * expfl
                              egbias = egbias + bias
                              dgdl = dgdl - ldelta*sigl2inv * bias
                              dgdfl = dgdfl
     &                                  - fldelta*sigf2inv * bias
                           end if
                        end do
                     end if
                     ihist = ostnext(ihist)
                  end do
               end if
            end do
         end if
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine egkernelinterpolate -- interpolate g kernel  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "egkernelinterpolate" evaluates the ost g kernel and its first
c     derivatives using bicubic Hermite interpolation on the grid
c
c
      subroutine egkernelinterpolate (egbias,dgdl,dgdfl)
      use dlmda
      use mutant
      use ost
      implicit none
      integer i,j
      integer ia,ja
      integer il0,if0
      real*8 egbias,dgdl,dgdfl
      real*8 flstart,flend
      real*8 l0,f0
      real*8 x,y
      real*8 x2,x3,y2,y3
      real*8 hxv(2),hxd(2),dhxv(2),dhxd(2)
      real*8 hyv(2),hyd(2),dhyv(2),dhyd(2)
      real*8 val,gl,gf,glf
c
c
c     require a complete interpolation cell in both dimensions
c
      egbias = 0.0d0
      dgdl = 0.0d0
      dgdfl = 0.0d0
      flstart = dble(1-fli0) * wflmda
      flend = dble(nflmda-fli0) * wflmda
      if (lambda .lt. 0.0d0 .or. lambda .gt. 1.0d0)
     &   return
      if (dedl .lt. flstart .or. dedl .gt. flend)  return
c
c     locate the lower-left grid point of the interpolation cell
c
      if (lambda .ge. 1.0d0) then
         il0 = nlmda - 1
      else
         il0 = int(lambda/wlmda) + 1
         il0 = max(1,min(il0,nlmda-1))
      end if
      if (dedl .ge. flend) then
         if0 = nflmda - 1
      else
         if0 = int((dedl-flstart)/wflmda) + 1
         if0 = max(1,min(if0,nflmda-1))
      end if
      l0 = dble(il0-1) * wlmda
      f0 = dble(if0-fli0) * wflmda
      x = (lambda-l0) / wlmda
      y = (dedl-f0) / wflmda
c
c     cubic Hermite basis functions and normalized derivatives
c
      x2 = x * x
      x3 = x2 * x
      y2 = y * y
      y3 = y2 * y
      hxv(1) = 2.0d0*x3 - 3.0d0*x2 + 1.0d0
      hxv(2) = -2.0d0*x3 + 3.0d0*x2
      hxd(1) = x3 - 2.0d0*x2 + x
      hxd(2) = x3 - x2
      dhxv(1) = 6.0d0*x2 - 6.0d0*x
      dhxv(2) = -6.0d0*x2 + 6.0d0*x
      dhxd(1) = 3.0d0*x2 - 4.0d0*x + 1.0d0
      dhxd(2) = 3.0d0*x2 - 2.0d0*x
      hyv(1) = 2.0d0*y3 - 3.0d0*y2 + 1.0d0
      hyv(2) = -2.0d0*y3 + 3.0d0*y2
      hyd(1) = y3 - 2.0d0*y2 + y
      hyd(2) = y3 - y2
      dhyv(1) = 6.0d0*y2 - 6.0d0*y
      dhyv(2) = -6.0d0*y2 + 6.0d0*y
      dhyd(1) = 3.0d0*y2 - 4.0d0*y + 1.0d0
      dhyd(2) = 3.0d0*y2 - 2.0d0*y
c
c     evaluate the bicubic Hermite patch and its derivatives
c
      do ia = 1, 2
         i = il0 + ia - 1
         do ja = 1, 2
            j = if0 + ja - 1
            val = gkernel(i,j)
            gl = wlmda * glkernel(i,j)
            gf = wflmda * gfkernel(i,j)
            glf = wlmda * wflmda * glfkernel(i,j)
            egbias = egbias
     &         + hxv(ia)*hyv(ja)*val
     &         + hxd(ia)*hyv(ja)*gl
     &         + hxv(ia)*hyd(ja)*gf
     &         + hxd(ia)*hyd(ja)*glf
            dgdl = dgdl
     &         + dhxv(ia)*hyv(ja)*val
     &         + dhxd(ia)*hyv(ja)*gl
     &         + dhxv(ia)*hyd(ja)*gf
     &         + dhxd(ia)*hyd(ja)*glf
            dgdfl = dgdfl
     &         + hxv(ia)*dhyv(ja)*val
     &         + hxd(ia)*dhyv(ja)*gl
     &         + hxv(ia)*dhyd(ja)*gf
     &         + hxd(ia)*dhyd(ja)*glf
         end do
      end do
      dgdl = dgdl / wlmda
      dgdfl = dgdfl / wflmda
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  function flambdabin -- get bin index for flambda  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "flambdabin" computes the bin index for dU/dlambda
c
c
      function flambdabin (dudl)
      use ost
      implicit none
      integer flambdabin
      real*8 dudl
c
c
c     set flambdabin value
c
      flambdabin = nint(dudl / wflmda) + fli0
      if (flambdabin .lt. 1)  flambdabin = 1
      if (flambdabin .gt. nflmda)  flambdabin = nflmda
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine resizeosthist -- resize ost histogram arrays  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "resizeosthist" doubles the storage for saved ost histogram
c     entries while preserving the previously saved data
c
c
      subroutine resizeosthist
      use dlmda
      use ost
      implicit none
      integer oldsize
      integer newsize
      integer, allocatable :: itmp(:)
      real*8, allocatable :: rtmp(:)
c
c
c     grow each array, keeping the saved entries and zeroing the rest
c
      oldsize = sizelmdahist
      newsize = 2 * oldsize
      allocate (itmp(newsize))
      itmp(1:oldsize) = osthist(1:oldsize)
      itmp(oldsize+1:newsize) = 0
      call move_alloc (itmp,osthist)
      allocate (itmp(newsize))
      itmp(1:oldsize) = lmdaihist(1:oldsize)
      itmp(oldsize+1:newsize) = 0
      call move_alloc (itmp,lmdaihist)
      allocate (itmp(newsize))
      itmp(1:oldsize) = ostnext(1:oldsize)
      itmp(oldsize+1:newsize) = 0
      call move_alloc (itmp,ostnext)
      allocate (rtmp(newsize))
      rtmp(1:oldsize) = lmdalhist(1:oldsize)
      rtmp(oldsize+1:newsize) = 0.0d0
      call move_alloc (rtmp,lmdalhist)
      allocate (rtmp(newsize))
      rtmp(1:oldsize) = lmdafhist(1:oldsize)
      rtmp(oldsize+1:newsize) = 0.0d0
      call move_alloc (rtmp,lmdafhist)
      allocate (rtmp(newsize))
      rtmp(1:oldsize) = osthhist(1:oldsize)
      rtmp(oldsize+1:newsize) = 0.0d0
      call move_alloc (rtmp,osthhist)
      allocate (rtmp(newsize))
      rtmp(1:oldsize) = ostwlhist(1:oldsize)
      rtmp(oldsize+1:newsize) = 0.0d0
      call move_alloc (rtmp,ostwlhist)
      allocate (rtmp(newsize))
      rtmp(1:oldsize) = ostwfhist(1:oldsize)
      rtmp(oldsize+1:newsize) = 0.0d0
      call move_alloc (rtmp,ostwfhist)
      sizelmdahist = newsize
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine buildostindex -- build ost histogram index  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "buildostindex" rebuilds the packed bin locations and
c     linked-list lookup table from the saved real centers
c
c
      subroutine buildostindex
      use dlmda
      use ost
      implicit none
      integer i,j
      integer ihist
      integer k
      integer ilmda,iflmda
      integer lmdabin,flambdabin
c
c
c     clear the bin heads and next links
c
      do i = 1, nlmda
         do j = 1, nflmda
            osthead(i,j) = 0
         end do
      end do
      do i = 1, sizelmdahist
         ostnext(i) = 0
      end do
c
c     insert each saved histogram entry into its lambda/flambda bin
c
      do ihist = 1, nlmdahist
         ilmda = lmdabin(lmdalhist(ihist))
         iflmda = flambdabin(lmdafhist(ihist))
         call ij_to_k (ilmda,iflmda,nlmda,k)
         osthist(ihist) = k
         ostnext(ihist) = osthead(ilmda,iflmda)
         osthead(ilmda,iflmda) = ihist
      end do
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine buildgkernel -- build the ost g kernel  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "buildgkernel" builds the ost g kernel by looping over
c     saved histogram sources and spreading their normalized
c     gaussian contribution to nearby target bins
c
c
      subroutine buildgkernel
      use dlmda
      use ost
      implicit none
      integer i,j
      integer ihist
c
c
c     zero out g kernel and its running maxima
c
      do i = 1, nlmda
         vkernelmax(i) = 0.0d0
         do j = 1, nflmda
            gkernel(i,j) = 0.0d0
         end do
      end do
c
c     loop over saved histogram sources
c
      do ihist = 1, nlmdahist
         call addkernelhist (ihist,.false.)
      end do
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine buildkernels -- build the ost kernels  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "buildkernels" builds both ost kernels by looping over saved
c     histogram sources and updating the f kernel accumulators
c
c
      subroutine buildkernels
      use dlmda
      use ost
      implicit none
      integer i,j
      integer ihist
c
c
c     zero out kernels and free energy accumulators
c
      do i = 1, nlmda
         lmdafmean(i) = 0.0d0
         lmdafsum(i) = 0.0d0
         lmdafwt(i) = 0.0d0
         vkernelmax(i) = 0.0d0
         do j = 1, nflmda
            gfkernel(i,j) = 0.0d0
            gkernel(i,j) = 0.0d0
            glfkernel(i,j) = 0.0d0
            glkernel(i,j) = 0.0d0
         end do
      end do
c
c     loop over saved histogram sources
c
      do ihist = 1, nlmdahist
         call addkernelhist (ihist,.true.)
      end do
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine updategkernel -- update the ost g kernel  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "updategkernel" updates the ost g kernel by spreading
c     the most recently saved histogram source to nearby target bins
c
c
      subroutine updategkernel
      use dlmda
      implicit none
c
c
c     add the most recently saved histogram source
c
      if (nlmdahist .gt. 0)  call addkernelhist (nlmdahist,.false.)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine updatekernels -- update the ost kernels  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "updatekernels" updates the ost kernels by spreading the most
c     recently saved histogram source and updating f accumulators
c
c
      subroutine updatekernels
      use dlmda
      implicit none
c
c
c     add the most recently saved histogram source
c
      if (nlmdahist .gt. 0)  call addkernelhist (nlmdahist,.true.)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine addkernelhist -- add one histogram source  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "addkernelhist" spreads one saved histogram source to the ost
c     g kernel using normalized gaussian biasing functions; when
c     "dofkern" is set it also updates the g kernel derivatives and
c     the f kernel accumulators, otherwise only the g kernel and its
c     running maxima
c
c
      subroutine addkernelhist (ihist,dofkern)
      use bath
      use dlmda
      use ost
      use math
      use units
      implicit none
      integer ihist
      integer k
      integer lsrc,fsrc
      integer ilmda,iflmda
      integer img,llog
      integer ilmda1,ilmda2
      integer iflmda1,iflmda2
      integer nlcut,nfcut
      integer nbincut
      real*8 sigl,sigf
      real*8 sigl2,sigf2
      real*8 pref,rt
      real*8 sourcel,sourcefl
      real*8 targetl,targetf
      real*8 ldelta,ldelta2
      real*8 fldelta,fldelta2
      real*8 e,expl,expfl
      logical dofkern
c
c
c     unpack saved source location and get gaussian parameters
c
      k = osthist(ihist)
      call k_to_ij (k,nlmda,lsrc,fsrc)
      sigl = ostwlhist(ihist)
      sigf = ostwfhist(ihist)
      sigl2 = sigl * sigl
      sigf2 = sigf * sigf
      pref = osthhist(ihist) / (2.0d0*pi*sigl*sigf)
      sourcefl = lmdafhist(ihist)
      rt = gasconst * kelvin
c
c     include target bins out to the requested gaussian cutoff
c
      nlcut = nbincut(oststdev*sigl,wlmda)
      nfcut = nbincut(oststdev*sigf,wflmda)
      iflmda1 = max(1,fsrc-nfcut)
      iflmda2 = min(nflmda,fsrc+nfcut)
c
c     use the physical lambda source plus its two mirror images
c
      do img = 1, 3
         if (img .eq. 1) then
            llog = lsrc
            sourcel = lmdalhist(ihist)
         else if (img .eq. 2) then
            llog = 2 - lsrc
            sourcel = -lmdalhist(ihist)
         else
            llog = 2*nlmda - lsrc
            sourcel = 2.0d0 - lmdalhist(ihist)
         end if
c
c     only target bins within the gaussian cutoff receive a contribution
c
         ilmda1 = max(1,llog-nlcut)
         ilmda2 = min(nlmda,llog+nlcut)
         do ilmda = ilmda1, ilmda2
            targetl = dble(ilmda-1) * wlmda
            ldelta = targetl - sourcel
            if (abs(ldelta) .gt. oststdev*sigl)  goto 10
            ldelta2 = ldelta*ldelta
            expl = exp(-ldelta2 / (2.0d0*sigl2))
            do iflmda = iflmda1, iflmda2
               targetf = dble(iflmda-fli0) * wflmda
               fldelta = targetf - sourcefl
               if (abs(fldelta) .gt. oststdev*sigf)  goto 20
               fldelta2 = fldelta*fldelta
               expfl = exp(-fldelta2 / (2.0d0*sigf2))
               e = pref * expl * expfl
               if (dofkern) then
                  call addkernelpoint (ilmda,iflmda,e,ldelta,
     &                                 fldelta,sigl2,sigf2,rt)
               else
                  gkernel(ilmda,iflmda) = gkernel(ilmda,iflmda) + e
                  vkernelmax(ilmda) = max(vkernelmax(ilmda),
     &                                    gkernel(ilmda,iflmda))
               end if
   20          continue
            end do
c
c     set the mean force of the row once all its points are added
c
            if (dofkern) then
               if (lmdafwt(ilmda) .eq. 0.0d0) then
                  lmdafmean(ilmda) = 0.0d0
               else
                  lmdafmean(ilmda) = lmdafsum(ilmda) / lmdafwt(ilmda)
               end if
            end if
   10       continue
         end do
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine addkernelpoint -- add one grid contribution  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "addkernelpoint" updates one g kernel cell and the associated
c     f kernel numerator and partition-function accumulator, given
c     the thermal energy "rt"; the caller sets the row mean force
c
c
      subroutine addkernelpoint (ilmda,iflmda,e,ldelta,
     &                           fldelta,sigl2,sigf2,rt)
      use dlmda
      use ost
      implicit none
      integer ilmda,iflmda
      real*8 e,rt
      real*8 ldelta,fldelta
      real*8 sigl2,sigf2
      real*8 flmda
      real*8 oldg,newg
      real*8 oldweight,newweight
      real*8 delweight
      real*8 scale,vmax
      real*8 dgdl,dgdfl,d2gdlfl
      logical raised
c
c
c     update the g kernel and adjust the f kernel accumulators
c
      oldg = gkernel(ilmda,iflmda)
      newg = oldg + e
c
c     the accumulators carry the largest bias in the row as a common
c     factor, so they are rescaled whenever this point raises it, and
c     the new point then has unit weight
c
      vmax = vkernelmax(ilmda)
      raised = (newg .gt. vmax)
      if (raised) then
         scale = exp((vmax-newg)/rt)
         lmdafsum(ilmda) = lmdafsum(ilmda) * scale
         lmdafwt(ilmda) = lmdafwt(ilmda) * scale
         vmax = newg
      end if
      if (oldg .eq. 0.0d0) then
         oldweight = 0.0d0
      else
         oldweight = exp((oldg-vmax)/rt)
      end if
      if (raised) then
         newweight = 1.0d0
      else
         newweight = exp((newg-vmax)/rt)
      end if
      delweight = newweight - oldweight
      flmda = dble(iflmda-fli0) * wflmda
      dgdl = -ldelta * e / sigl2
      dgdfl = -fldelta * e / sigf2
      d2gdlfl = ldelta * fldelta * e / (sigl2*sigf2)
      gfkernel(ilmda,iflmda) = gfkernel(ilmda,iflmda) + dgdfl
      gkernel(ilmda,iflmda) = newg
      vkernelmax(ilmda) = max(vkernelmax(ilmda),newg)
      glfkernel(ilmda,iflmda) = glfkernel(ilmda,iflmda) + d2gdlfl
      glkernel(ilmda,iflmda) = glkernel(ilmda,iflmda) + dgdl
      lmdafsum(ilmda) = lmdafsum(ilmda) + flmda*delweight
      lmdafwt(ilmda) = lmdafwt(ilmda) + delweight
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine buildfkernel -- build the ost f kernel  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "buildfkernel" builds the ost f kernel by computing the
c     biased ensemble average of dE/dlambda at each lambda bin
c
c
      subroutine buildfkernel
      use bath
      use dlmda
      use ost
      use units
      implicit none
      integer ilmda,iflmda
      real*8 avgflambda
      real*8 partfunc
      real*8 flmda
      real*8 vmax
      real*8 weight
c
c
c     loop over lambda bins
c
      do ilmda = 1, nlmda
         avgflambda = 0.0d0
         partfunc = 0.0d0
c
c     the largest bias in the row is factored out of both sums, where
c     it cancels exactly, to keep the exponential from overflowing
c
         vmax = vkernelmax(ilmda)
c
c     loop over flambda bins with nonzero g kernel support
c
         do iflmda = 1, nflmda
            if (gkernel(ilmda,iflmda) .ne. 0.0d0) then
               flmda = dble(iflmda-fli0) * wflmda
               weight = exp((gkernel(ilmda,iflmda)-vmax)
     &                         /(gasconst*kelvin))
               avgflambda = avgflambda + flmda*weight
               partfunc = partfunc + weight
            end if
         end do
c
c     set mean force for this lambda bin
c
         if (partfunc .eq. 0.0d0) then
            lmdafmean(ilmda) = 0.0d0
         else
            lmdafmean(ilmda) = avgflambda / partfunc
         end if
      end do
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine initostfile -- open the ost history file  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "initostfile" opens the external file that will hold the ost
c     bias history and writes its initial header; a new version is
c     claimed so a rerun cannot append onto an earlier history
c
c
      subroutine initostfile
      use dlmda
      use files
      use iounit
      implicit none
      integer ihis
      integer freeunit
      integer trimtext
c
c
c     return if ost has not been initialized
c
      if (.not. use_ost)  return
c
c     open a new file, keeping any output from a previous run
c
      ihis = freeunit ()
      lmdasavefile = filename(1:leng)//'.ost'
      call version (lmdasavefile,'new')
      open (unit=ihis,file=lmdasavefile,status='new',
     &      access='stream',form='formatted')
      call prtost (ihis)
      close (unit=ihis)
      nlmdasave = nlmdahist
c
c     report the name of the file holding the bias history
c
      write (iout,10)  lmdasavefile(1:trimtext(lmdasavefile))
   10 format (/,' OST  --  Bias History Written To  ',a)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine saveost -- output ost history information  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "saveost" appends to the external .ost file the bias history
c     accumulated since the previous call, so a running simulation
c     can be monitored; the fixed-size header is updated in place
c
c
      subroutine saveost
      use dlmda
      use ost
      implicit none
      integer ihis
      integer freeunit
c
c
c     return if ost has not been initialized
c
      if (.not. use_ost)  return
      if (.not. allocated(osthist))  return
c
c     append the history recorded since the previous call before
c     updating the header, so an interrupted write leaves the
c     previous history count as the last complete checkpoint
c
      ihis = freeunit ()
      if (nlmdahist .gt. nlmdasave) then
         open (unit=ihis,file=lmdasavefile,status='old',
     &         access='stream',form='formatted',position='append')
         call prtosthist (ihis,nlmdasave+1,nlmdahist)
         close (unit=ihis)
      end if
      open (unit=ihis,file=lmdasavefile,status='old',
     &      access='stream',form='unformatted',position='rewind')
      call updosthead (ihis)
      close (unit=ihis)
      nlmdasave = nlmdahist
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine rdost -- input ost history information  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "rdost" reads a saved ost bias history from the external
c     .ost file, then rebuilds the histogram lookup and bias
c     kernels; used for analysis, not to restart a simulation
c
c
      subroutine rdost
      use dlmda
      use iounit
      use ost
      implicit none
      integer ihis
      integer freeunit
      integer trimtext
      logical exist
      character*240 ostfile
      character*240 title
c
c
c     return unless an ost history file is present
c
      ostfile = lmdasavefile
      inquire (file=ostfile,exist=exist)
      if (.not. exist)  return
      ihis = freeunit ()
      open (unit=ihis,file=ostfile,status='old')
      rewind (unit=ihis)
c
c     read the header, which must be that of an ost history
c
      call rdbiashead (ihis,ostfile,title)
      if (index(title,osttitle(1:trimtext(osttitle))) .eq. 0) then
         close (unit=ihis)
         write (iout,10)  ostfile(1:trimtext(ostfile))
   10    format (/,' RDOST  --  File is not an OST History File',
     &           /,'           File Name :  ',a)
         call fatal
      end if
c
c     read the gaussian history entries and rebuild the kernels
c
      call rdosthist (ihis,ostfile)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine rdosthist -- input ost history entries  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "rdosthist" reads the gaussian history entries that follow the
c     header of an open .ost file, then closes the file and rebuilds
c     the histogram lookup and bias kernels
c
c
      subroutine rdosthist (ihis,ostfile)
      use dlmda
      use inform
      use iounit
      use ost
      implicit none
      integer i,ihis
      integer ihist
      integer osthist0
      integer trimtext
      real*8 efreetot
      character*240 record
      character*(*) ostfile
c
c
c     reallocate ost arrays to match the history file
c
      if (allocated(osthhist))  deallocate (osthhist)
      if (allocated(osthist))  deallocate (osthist)
      if (allocated(lmdaihist))  deallocate (lmdaihist)
      if (allocated(osthead))  deallocate (osthead)
      if (allocated(ostnext))  deallocate (ostnext)
      if (allocated(lmdallist))  deallocate (lmdallist)
      if (allocated(lmdaflist))  deallocate (lmdaflist)
      if (allocated(lmdalhist))  deallocate (lmdalhist)
      if (allocated(lmdafhist))  deallocate (lmdafhist)
      if (allocated(ostwlhist))  deallocate (ostwlhist)
      if (allocated(ostwfhist))  deallocate (ostwfhist)
      if (allocated(lmdafmean))  deallocate (lmdafmean)
      if (allocated(lmdafsum))  deallocate (lmdafsum)
      if (allocated(gfkernel))  deallocate (gfkernel)
      if (allocated(gkernel))  deallocate (gkernel)
      if (allocated(glfkernel))  deallocate (glfkernel)
      if (allocated(glkernel))  deallocate (glkernel)
      if (allocated(lmdafwt))  deallocate (lmdafwt)
      if (allocated(vkernelmax))  deallocate (vkernelmax)
      allocate (osthhist(sizelmdahist))
      allocate (osthist(sizelmdahist))
      allocate (lmdaihist(sizelmdahist))
      allocate (osthead(nlmda,nflmda))
      allocate (ostnext(sizelmdahist))
      allocate (lmdallist(lmdaintv))
      allocate (lmdaflist(lmdaintv))
      allocate (lmdalhist(sizelmdahist))
      allocate (lmdafhist(sizelmdahist))
      allocate (ostwlhist(sizelmdahist))
      allocate (ostwfhist(sizelmdahist))
      allocate (lmdafmean(nlmda))
      allocate (lmdafsum(nlmda))
      allocate (gfkernel(nlmda,nflmda))
      allocate (gkernel(nlmda,nflmda))
      allocate (glfkernel(nlmda,nflmda))
      allocate (glkernel(nlmda,nflmda))
      allocate (lmdafwt(nlmda))
      allocate (vkernelmax(nlmda))
c
c     initialize ost arrays
c
      do i = 1, nlmda
         lmdafmean(i) = 0.0d0
         lmdafsum(i) = 0.0d0
         lmdafwt(i) = 0.0d0
         vkernelmax(i) = 0.0d0
         do ihist = 1, nflmda
            gfkernel(i,ihist) = 0.0d0
            gkernel(i,ihist) = 0.0d0
            glfkernel(i,ihist) = 0.0d0
            glkernel(i,ihist) = 0.0d0
            osthead(i,ihist) = 0
         end do
      end do
      do i = 1, lmdaintv
         lmdallist(i) = 0.0d0
         lmdaflist(i) = 0.0d0
      end do
      do i = 1, sizelmdahist
         osthist(i) = 0
         lmdaihist(i) = 0
         ostnext(i) = 0
         lmdalhist(i) = 0.0d0
         lmdafhist(i) = 0.0d0
         osthhist(i) = 0.0d0
         ostwlhist(i) = 0.0d0
         ostwfhist(i) = 0.0d0
      end do
c
c     start the gaussian width maxima from the widths saved with
c     each gaussian
c
      maxwlhist = 0.0d0
      maxwfhist = 0.0d0
c
c     read saved gaussian histogram entries
c
      do ihist = 1, nlmdahist
         read (ihis,10,err=90,end=90)  record
         read (record,*,err=90,end=90)  lmdaihist(ihist),osthist0,
     &      lmdalhist(ihist),lmdafhist(ihist),osthhist(ihist),
     &      ostwlhist(ihist),ostwfhist(ihist)
         osthist(ihist) = osthist0
         maxwlhist = max(maxwlhist,ostwlhist(ihist))
         maxwfhist = max(maxwfhist,ostwfhist(ihist))
      end do
   10 format (a240)
      close (unit=ihis)
c
c     rebuild lookup table and kernels from the saved gaussians
c
      call buildostindex
      if (fastkernel) then
         call buildkernels
      else
         call buildgkernel
         call buildfkernel
      end if
      lmdadeltag = efreetot()
      nlmdasave = nlmdahist
      if (debug) then
         write (iout,20)  ostfile(1:trimtext(ostfile))
   20    format (/,' Reading OST Bias from :  ',a)
      end if
      return
c
c     malformed history file
c
   90 continue
      close (unit=ihis)
      write (iout,30)  ostfile(1:trimtext(ostfile))
   30 format (/,' RDOST  --  Error while Reading OST History File',
     &        /,'            File Name :  ',a)
      call fatal
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine prtost -- output of ost energy and histogram  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prtost" writes current ost history state
c
c
      subroutine prtost (ihis)
      use dlmda
      implicit none
      integer ihis
c
c     write the header followed by all current history entries
c
      call prtosthead (ihis)
      call prtosthist (ihis,1,nlmdahist)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine prtosthead  --  output ost history header  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "prtosthead" writes the fixed-size current ost history header
c
c
      subroutine prtosthead (ihis)
      use ost
      implicit none
      integer ihis
c
c
c     write the header under the ost title
c
      call prtbiashead (ihis,osttitle,ostlabel)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine updosthead  --  update ost history header  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "updosthead" overwrites the fixed-size ost header in place;
c     unformatted stream output avoids truncating the appended history
c
c
      subroutine updosthead (ihis)
      use ost
      implicit none
      integer ihis
c
c
c     update the header under the ost title
c
      call updbiashead (ihis,osttitle,ostlabel)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine prtosthist  --  output ost history range  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "prtosthist" writes a requested range of gaussian histories
c
c
      subroutine prtosthist (ihis,ifirst,ilast)
      use dlmda
      use ost
      implicit none
      integer ihis
      integer ifirst
      integer ilast
      integer ihist
c
c     write saved gaussian histogram entries
c
      do ihist = ifirst, ilast
         write (ihis,10)  lmdaihist(ihist),osthist(ihist),
     &      lmdalhist(ihist),lmdafhist(ihist),osthhist(ihist),
     &      ostwlhist(ihist),ostwfhist(ihist)
      end do
   10 format (2i12,5d26.16)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine initmetafile  --  open metadyn history file  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "initmetafile" opens the external file that will hold the
c     metadynamics bias history and writes its initial header; a new
c     version is claimed so a rerun cannot append onto an earlier one
c
c
      subroutine initmetafile
      use dlmda
      use files
      use iounit
      use ost
      implicit none
      integer ihis
      integer freeunit
      integer trimtext
c
c
c     return if metadynamics has not been initialized
c
      if (.not. use_meta)  return
c
c     open a new file, keeping any output from a previous run
c
      ihis = freeunit ()
      metasavefile = filename(1:leng)//'.meta'
      call version (metasavefile,'new')
      open (unit=ihis,file=metasavefile,status='new',
     &      access='stream',form='formatted')
      call prtmeta (ihis)
      close (unit=ihis)
      nmethistsave = nmetahist
c
c     report the name of the file holding the bias history
c
      write (iout,10)  metasavefile(1:trimtext(metasavefile))
   10 format (/,' METADYN  --  Bias History Written To  ',a)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine savemeta -- output metadynamics history  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "savemeta" appends to the external .meta file the bias history
c     accumulated since the previous call, so a running simulation
c     can be monitored; the fixed-size header is updated in place
c
c
      subroutine savemeta
      use dlmda
      use ost
      implicit none
      integer ihis
      integer freeunit
c
c
c     return if metadynamics has not been initialized
c
      if (.not. use_meta)  return
      if (.not. allocated(metalhist))  return
c
c     append the history recorded since the previous call before
c     updating the header, so an interrupted write leaves the
c     previous history count as the last complete checkpoint
c
      ihis = freeunit ()
      if (nmetahist .gt. nmethistsave) then
         open (unit=ihis,file=metasavefile,status='old',
     &         access='stream',form='formatted',position='append')
         call prtmetahist (ihis,nmethistsave+1,nmetahist)
         close (unit=ihis)
      end if
      open (unit=ihis,file=metasavefile,status='old',
     &      access='stream',form='unformatted',position='rewind')
      call updmetahead (ihis)
      close (unit=ihis)
      nmethistsave = nmetahist
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine rdmetahist  --  input metadynamics history  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "rdmetahist" reads the gaussian history entries that follow the
c     header of an open .meta file, then closes the file and rebuilds
c     the bin lookup and the metadynamics bias grid
c
c
      subroutine rdmetahist (ihis,metafile)
      use dlmda
      use inform
      use iounit
      use ost
      implicit none
      integer i,ihis
      integer ihist
      integer trimtext
      real*8 metadeltag
      character*240 record
      character*(*) metafile
c
c
c     reallocate metadynamics arrays to match the history file
c
      if (allocated(metalhist))  deallocate (metalhist)
      if (allocated(metahhist))  deallocate (metahhist)
      if (allocated(metawhist))  deallocate (metawhist)
      if (allocated(metaihist))  deallocate (metaihist)
      if (allocated(metanext))  deallocate (metanext)
      if (allocated(metahead))  deallocate (metahead)
      if (allocated(lmdallist))  deallocate (lmdallist)
      if (allocated(vmetagrid))  deallocate (vmetagrid)
      if (allocated(dvmetagrid))  deallocate (dvmetagrid)
      allocate (metalhist(sizemetahist))
      allocate (metahhist(sizemetahist))
      allocate (metawhist(sizemetahist))
      allocate (metaihist(sizemetahist))
      allocate (metanext(sizemetahist))
      allocate (metahead(nlmda))
      allocate (lmdallist(lmdaintv))
      allocate (vmetagrid(nlmda))
      allocate (dvmetagrid(nlmda))
c
c     initialize metadynamics arrays
c
      do i = 1, nlmda
         vmetagrid(i) = 0.0d0
         dvmetagrid(i) = 0.0d0
      end do
      do i = 1, lmdaintv
         lmdallist(i) = 0.0d0
      end do
      do i = 1, sizemetahist
         metalhist(i) = 0.0d0
         metahhist(i) = 0.0d0
         metawhist(i) = 0.0d0
         metaihist(i) = 0
      end do
c
c     read saved metadynamics gaussian entries
c
      do ihist = 1, nmetahist
         read (ihis,10,err=90,end=90)  record
         read (record,*,err=90,end=90)  metaihist(ihist),
     &      metalhist(ihist),metahhist(ihist),metawhist(ihist)
      end do
   10 format (a240)
      close (unit=ihis)
c
c     take the gaussian width maximum from the widths saved with each
c     gaussian, then rebuild the bin lookup and the grid from them
c
      maxwlhist = 0.0d0
      do ihist = 1, nmetahist
         maxwlhist = max(maxwlhist,metawhist(ihist))
      end do
      call buildmetaindex
      do ihist = 1, nmetahist
         call addmetagrid (ihist)
      end do
      lmdadeltag = metadeltag()
      nmethistsave = nmetahist
      if (debug) then
         write (iout,20)  metafile(1:trimtext(metafile))
   20    format (/,' Reading Metadynamics Bias from :  ',a)
      end if
      return
c
c     malformed history file
c
   90 continue
      close (unit=ihis)
      write (iout,30)  metafile(1:trimtext(metafile))
   30 format (/,' RDMETAHIST  --  Error while Reading Metadynamics',
     &           ' History File',
     &        /,'                 File Name :  ',a)
      call fatal
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine prtmeta -- output metadynamics history  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "prtmeta" writes current metadynamics history state
c
c
      subroutine prtmeta (ihis)
      use ost
      implicit none
      integer ihis
c
c
c     write the header followed by all current history entries
c
      call prtmetahead (ihis)
      call prtmetahist (ihis,1,nmetahist)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine prtmetahead  --  output metadynamics header  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "prtmetahead" writes the fixed-size metadynamics history header
c
c
      subroutine prtmetahead (ihis)
      use ost
      implicit none
      integer ihis
c
c
c     write the header under the metadynamics title
c
      call prtbiashead (ihis,metatitle,metalabel)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine updmetahead  --  update metadynamics header  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "updmetahead" overwrites the fixed-size metadynamics header in
c     place; unformatted stream output avoids truncating the history
c
c
      subroutine updmetahead (ihis)
      use ost
      implicit none
      integer ihis
c
c
c     update the header under the metadynamics title
c
      call updbiashead (ihis,metatitle,metalabel)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine prtmetahist  --  output metadynamics history  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prtmetahist" writes a requested range of metadynamics gaussians
c
c
      subroutine prtmetahist (ihis,ifirst,ilast)
      use ost
      implicit none
      integer ihis
      integer ifirst
      integer ilast
      integer ihist
c
c     write saved metadynamics gaussian entries
c
      do ihist = ifirst, ilast
         write (ihis,10)  metaihist(ihist),metalhist(ihist),
     &      metahhist(ihist),metawhist(ihist)
      end do
   10 format (i12,3d26.16)
      return
      end
c
c
c     ####################################################
c     ##                                                ##
c     ##  function nbincut -- bins spanned by a cutoff  ##
c     ##                                                ##
c     ####################################################
c
c
c     "nbincut" returns the smallest number of bins of the given
c     width that together span the given cutoff distance
c
c
      function nbincut (cut,width)
      implicit none
      integer nbincut
      real*8 cut,width
c
c
c     round the cutoff in bin widths up to a whole number
c
      nbincut = int(cut/width)
      if (dble(nbincut)*width .lt. cut)  nbincut = nbincut + 1
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine lmdaimgbin -- lambda bin of an image  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "lmdaimgbin" maps a lambda bin index "lcenter", which may lie
c     beyond either end of the lambda range, onto the bin "lcount"
c     whose gaussians have their mirror image in it, and returns in
c     "nimg" the number of images to count there, two in each end
c     bin for the gaussian and its reflection, otherwise one
c
c
      subroutine lmdaimgbin (lcenter,lcount,nimg)
      use dlmda
      implicit none
      integer lcenter,lcount
      integer nimg
c
c
c     reflect the bin index at the lambda end points
c
      lcount = lcenter
      if (lcount .lt. 1) then
         lcount = 2 - lcount
      else if (lcount .gt. nlmda) then
         lcount = 2*nlmda - lcount
      end if
      nimg = 1
      if (lcenter.eq.1 .or. lcenter.eq.nlmda)  nimg = 2
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  function lmdaimgsrc -- lambda center of an image  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "lmdaimgsrc" returns the lambda center of image "img" of a
c     gaussian centered at "lsrc", as counted in lambda bin "lcenter"
c     by "lmdaimgbin"; beyond the ends it is the reflected center,
c     and in an end bin the second image is the reflected one
c
c
      function lmdaimgsrc (lcenter,img,lsrc)
      use dlmda
      implicit none
      integer lcenter,img
      real*8 lmdaimgsrc
      real*8 lsrc
c
c
c     choose the physical center or its reflection at an end point
c
      if (lcenter .lt. 1) then
         lmdaimgsrc = -lsrc
      else if (lcenter .gt. nlmda) then
         lmdaimgsrc = 2.0d0 - lsrc
      else if (img.eq.2 .and. lcenter.eq.1) then
         lmdaimgsrc = -lsrc
      else if (img.eq.2 .and. lcenter.eq.nlmda) then
         lmdaimgsrc = 2.0d0 - lsrc
      else
         lmdaimgsrc = lsrc
      end if
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine ij_to_k -- convert 2D index to 1D index  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "ij_to_k" converts the 2D indices i and j into the
c     corresponding 1D index k
c
c
      subroutine ij_to_k (i,j,nrow,k)
      implicit none
      integer i,j
      integer nrow
      integer k
c
c
c     convert 2D indices into a packed 1D index
c
      k = i + (j-1)*nrow
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine k_to_ij -- convert 1D index to 2D index  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "k_to_ij" converts the 1D index k into the corresponding
c     2D indices i and j
c
c
      subroutine k_to_ij (k,nrow,i,j)
      implicit none
      integer k
      integer nrow
      integer i,j
c
c
c     convert a packed 1D index into 2D indices
c
      i = mod(k-1,nrow) + 1
      j = (k-1)/nrow + 1
      return
      end
