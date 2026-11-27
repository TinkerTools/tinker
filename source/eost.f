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
c     evaluate the ost g kernel bias and free energy f term, where
c     the f term uses the mean force with unsampled bins filled
c
      if (ostinterpol) then
         call egkernelinterpolate (egbias,dgdl,dgdfl)
      else
         call egkernel (egbias,dgdl,dgdfl)
      end if
      call efreebias (lambda,eostlmda,dfdl)
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
c     ############################################################
c     ##                                                        ##
c     ##  subroutine eostvbias -- ost bias energy from history  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "eostvbias" finds the orthogonal space tempering bias energy
c     at the current lambda and unbiased dU/dlambda from the history
c     as it now stands, the same energy that "eostbias" adds, so a
c     deposit made after the forces of a step were found can be taken
c     into the energy of that step without a new gradient; no saved
c     value is changed
c
c
      subroutine eostvbias (vbias)
      use mutant
      use ost
      implicit none
      real*8 vbias
      real*8 egbias,dgdl,dgdfl
      real*8 eostlmda,dfdl
c
c
c     the g kernel bias less the free energy term, as in eostbias
c
      if (ostinterpol) then
         call egkernelinterpolate (egbias,dgdl,dgdfl)
      else
         call egkernel (egbias,dgdl,dgdfl)
      end if
      call efreebias (lambda,eostlmda,dfdl)
      vbias = egbias - eostlmda
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
c
c     add the new gaussian to the sample mean force, weighted by the
c     bias it was sampled under, before it joins the kernels
c
      call addflatsample (nlmdahist)
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
c     an untempered run always deposits at the full height; the
c     tempering scale is the bath temperature, not that of lambda
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
c     by direct summation over the whole saved history
c
c
      subroutine egkernel (egbias,dgdl,dgdfl)
      use dlmda
      use mutant
      implicit none
      real*8 egbias,dgdl,dgdfl
c
c
      call egkernelpt (lambda,dedl,nlmdahist,egbias,dgdl,dgdfl)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine egkernelpt -- ost g energy at a given point  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "egkernelpt" computes the ost biasing g kernel energy and
c     derivatives at lambda "lmda" and dU/dlambda "dudl" by direct
c     summation over the saved gaussians with history index up to
c     "maxhist", so the bias as it stood before a later gaussian was
c     added can be found from a history whose lookup holds them all
c
c
      subroutine egkernelpt (lmda,dudl,maxhist,egbias,dgdl,dgdfl)
      use dlmda
      use math
      use ost
      implicit none
      integer maxhist
      integer ihist
      integer ilmda,iflmda
      integer klmda,kflmda
      integer lcenter,lcount
      integer flcenter
      integer nlcut,nfcut
      integer img,nimg
      integer lmdabin
      integer nbincut
      real*8 lmda,dudl
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
c     get the lambda and flambda bin indices of the point
c
      ilmda = lmdabin(lmda)
      iflmda = nint(dudl / wflmda) + fli0
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
                     if (ihist .gt. maxhist)  goto 10
                     sigl = ostwlhist(ihist)
                     sigf = ostwfhist(ihist)
                     sourcefl = lmdafhist(ihist)
                     fldelta = dudl - sourcefl
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
                           ldelta = lmda - sourcel
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
   10                continue
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
c     derivatives using bicubic Hermite interpolation on the grid at
c     the current lambda and unbiased dU/dlambda values
c
c
      subroutine egkernelinterpolate (egbias,dgdl,dgdfl)
      use dlmda
      use mutant
      implicit none
      real*8 egbias,dgdl,dgdfl
c
c
      call egkernelinterppt (lambda,dedl,egbias,dgdl,dgdfl)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine egkernelinterppt -- interpolate g at a point  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "egkernelinterppt" evaluates the ost g kernel and its first
c     derivatives at lambda "lmda" and dU/dlambda "dudl" using bicubic
c     Hermite interpolation on the grid as it currently stands
c
c
      subroutine egkernelinterppt (lmda,dudl,egbias,dgdl,dgdfl)
      use dlmda
      use ost
      implicit none
      integer i,j
      integer ia,ja
      integer il0,if0
      real*8 lmda,dudl
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
      if (lmda .lt. 0.0d0 .or. lmda .gt. 1.0d0)
     &   return
      if (dudl .lt. flstart .or. dudl .gt. flend)  return
c
c     locate the lower-left grid point of the interpolation cell
c
      if (lmda .ge. 1.0d0) then
         il0 = nlmda - 1
      else
         il0 = int(lmda/wlmda) + 1
         il0 = max(1,min(il0,nlmda-1))
      end if
      if (dudl .ge. flend) then
         if0 = nflmda - 1
      else
         if0 = int((dudl-flstart)/wflmda) + 1
         if0 = max(1,min(if0,nflmda-1))
      end if
      l0 = dble(il0-1) * wlmda
      f0 = dble(if0-fli0) * wflmda
      x = (lmda-l0) / wlmda
      y = (dudl-f0) / wflmda
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
c     "buildgkernel" builds the ost g kernel and its derivatives by
c     looping over saved histogram sources and spreading their
c     normalized gaussian contribution to nearby target bins, along
c     with the sample mean force when OST-FLATTEN selects it, while
c     "buildfkernel" then gives the kernel mean force
c
c
      subroutine buildgkernel
      use dlmda
      use ost
      implicit none
      integer ihist
c
c
c     zero out g kernel, its derivatives and running maxima, and the
c     mean force of the lambda bins
c
      call clearkernels
c
c     loop over saved histogram sources, where each one joins the
c     sample mean force before it joins the g kernel
c
      do ihist = 1, nlmdahist
         call addflatsample (ihist)
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
c     "buildkernels" builds the ost g kernel, its derivatives and the
c     mean force OST-FLATTEN selects, kernel or sample, by looping over
c     the saved histogram sources with the fused kernel updates
c
c
      subroutine buildkernels
      use dlmda
      use ost
      implicit none
      integer ihist
c
c
c     zero out kernels and free energy accumulators
c
      call clearkernels
c
c     loop over saved histogram sources, where each one joins the
c     sample mean force before it joins the kernels
c
      do ihist = 1, nlmdahist
         call addflatsample (ihist)
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
c     g kernel, its derivatives and running maxima using normalized
c     gaussian biasing functions; when "dofkern" is set it also
c     updates the f kernel accumulators, unless OST-FLATTEN selects
c     the sample mean force
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
      integer iflo,ifhi
      integer nlcut,nfcut
      integer nbincut
      real*8 sigl,sigf
      real*8 sigl2,sigf2
      real*8 pref,rt,rtinv
      real*8 fcut,sigf2inv
      real*8 lfac,ffac
      real*8 sourcel,sourcefl
      real*8 targetl,targetf
      real*8 ldelta,ldelta2
      real*8 fldelta,fldelta2
      real*8 e,expl
      real*8 flmda
      real*8 oldg,newg
      real*8 oldweight,newweight
      real*8 delweight
      real*8 scale,vmax,vtop
      real*8 dgdl,dgdfl,d2gdlfl
      real*8 fexp(nflmda)
      logical dofkern
      logical dofmean
      logical raised
c
c
c     the f kernel accumulators hold the kernel mean force, so they
c     are kept only when that mean force is selected
c
      dofmean = (dofkern .and. .not.ostflatsamp)
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
c
c     the f kernel weights undo the bias on the atoms, so they use
c     the bath temperature, never the lambda particle temperature
c
      rt = gasconst * kelvin
      rtinv = 1.0d0 / rt
      sigf2inv = 1.0d0 / sigf2
c
c     include target bins out to the requested gaussian cutoff
c
      nlcut = nbincut(oststdev*sigl,wlmda)
      nfcut = nbincut(oststdev*sigf,wflmda)
      iflmda1 = max(1,fsrc-nfcut)
      iflmda2 = min(nflmda,fsrc+nfcut)
c
c     the gaussian factors along flambda are the same for every
c     lambda target, so get them once for all lambda bins; since
c     the flambda distance grows with the bin index, the bins
c     within the gaussian cutoff form one range
c
      fcut = oststdev * sigf
      iflo = iflmda2 + 1
      ifhi = iflmda1 - 1
      do iflmda = iflmda1, iflmda2
         targetf = dble(iflmda-fli0) * wflmda
         fldelta = targetf - sourcefl
         if (abs(fldelta) .le. fcut) then
            iflo = min(iflo,iflmda)
            ifhi = iflmda
            fldelta2 = fldelta*fldelta
            fexp(iflmda) = exp(-fldelta2 / (2.0d0*sigf2))
         end if
      end do
      iflmda1 = iflo
      iflmda2 = ifhi
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
            lfac = -ldelta / sigl2
            do iflmda = iflmda1, iflmda2
               e = pref * expl * fexp(iflmda)
               oldg = gkernel(ilmda,iflmda)
               newg = oldg + e
               vmax = vkernelmax(ilmda)
               raised = (newg .gt. vmax)
               vtop = max(vmax,newg)
               flmda = dble(iflmda-fli0) * wflmda
c
c     update the f kernel accumulators, which carry the largest bias
c     in the row as a common factor, so they are rescaled whenever
c     this point raises it, and the new point then has unit weight;
c     a cell with no bias carries no weight
c
               if (dofmean) then
                  if (raised) then
                     scale = exp((vmax-newg)*rtinv)
                     lmdafsum(ilmda) = lmdafsum(ilmda) * scale
                     lmdafwt(ilmda) = lmdafwt(ilmda) * scale
                  end if
                  if (oldg .eq. 0.0d0) then
                     oldweight = 0.0d0
                  else
                     oldweight = exp((oldg-vtop)*rtinv)
                  end if
                  if (newg .eq. 0.0d0) then
                     newweight = 0.0d0
                  else if (raised) then
                     newweight = 1.0d0
                  else
                     newweight = exp((newg-vtop)*rtinv)
                  end if
                  delweight = newweight - oldweight
                  lmdafsum(ilmda) = lmdafsum(ilmda) + flmda*delweight
                  lmdafwt(ilmda) = lmdafwt(ilmda) + delweight
               end if
c
c     update the g kernel, its running maximum and the derivatives
c     that the interpolation reads
c
               ffac = (sourcefl-flmda) * sigf2inv
               dgdl = lfac * e
               dgdfl = ffac * e
               d2gdlfl = lfac * dgdfl
               gkernel(ilmda,iflmda) = newg
               vkernelmax(ilmda) = vtop
               gfkernel(ilmda,iflmda) = gfkernel(ilmda,iflmda) + dgdfl
               glkernel(ilmda,iflmda) = glkernel(ilmda,iflmda) + dgdl
               glfkernel(ilmda,iflmda) = glfkernel(ilmda,iflmda)
     &                                      + d2gdlfl
            end do
c
c     set the mean force of the row once all its points are added
c
            if (dofmean) then
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
c     #########################################################
c     ##                                                     ##
c     ##  subroutine buildfkernel -- build the ost f kernel  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "buildfkernel" builds the ost f kernel by computing the
c     biased ensemble average of dE/dlambda at each lambda bin, and
c     keeps its weighted sum and total weight in the same scaled
c     form as "addkernelhist", so a row has data exactly when its
c     total weight is positive; it does nothing when OST-FLATTEN
c     selects the sample mean force
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
c     the f kernel holds the kernel mean force alone
c
      if (ostflatsamp)  return
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
c     loop over flambda bins with nonzero g kernel support; the
c     weights undo the bias on the atoms at the bath temperature
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
c     keep the accumulators as the fused kernel build does, and set
c     the mean force for this lambda bin
c
         lmdafsum(ilmda) = avgflambda
         lmdafwt(ilmda) = partfunc
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
c     ############################################################
c     ##                                                        ##
c     ##  function ostgdep -- ost bias before a saved gaussian  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "ostgdep" finds the ost g bias at the lambda and dU/dlambda
c     center of saved gaussian "ihist" from the gaussians saved before
c     it, the bias under which that interval was sampled; it uses the
c     same evaluation as the dynamics, bicubic interpolation of the
c     grid, which must not yet hold gaussian "ihist", or the direct
c     sum over the earlier gaussians
c
      function ostgdep (ihist)
      use dlmda
      use ost
      implicit none
      integer ihist
      real*8 ostgdep
      real*8 g,dgdl,dgdfl
c
c
c     interpolate the grid, or sum the earlier gaussians directly
c
      if (ostinterpol) then
         call egkernelinterppt (lmdalhist(ihist),lmdafhist(ihist),
     &                          g,dgdl,dgdfl)
      else
         call egkernelpt (lmdalhist(ihist),lmdafhist(ihist),ihist-1,
     &                    g,dgdl,dgdfl)
      end if
      ostgdep = g
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine addflatsample -- add one sample to mean force  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "addflatsample" adds saved gaussian "ihist" to the sample mean
c     force, held in the mean force of the lambda bins, when
c     OST-FLATTEN selects it, and otherwise does nothing; its
c     dU/dlambda center counts in the lambda bin of
c     its lambda center with weight exp(g/RT), and g is the bias the
c     dynamics felt there before this gaussian was added; the weights
c     undo that bias on the atoms, so they use the bath temperature
c
c     note the sums of a bin are kept relative to the largest g they
c     hold and rescaled when a sample raises it, so no weight exceeds
c     one and the exponential cannot overflow; the sums are those of
c     the mean force of the lambda bins, which belong to one choice of
c     OST-FLATTEN at a time, as explained in "clearkernels"
c
c
      subroutine addflatsample (ihist)
      use bath
      use dlmda
      use ost
      use units
      implicit none
      integer ihist
      integer ilmda
      integer lmdabin
      real*8 g,w,rt
      real*8 scale
      real*8 ostgdep
c
c
c     get the lambda bin and the bias at the gaussian center
c
      if (.not. ostflatsamp)  return
      rt = gasconst * kelvin
      ilmda = lmdabin(lmdalhist(ihist))
      g = ostgdep(ihist)
c
c     the first sample of a bin sets its reference, and a larger g
c     becomes the reference after the old sums are rescaled
c
      w = 1.0d0
      if (lmdafwt(ilmda) .eq. 0.0d0) then
         ostflatgmax(ilmda) = g
      else if (g .gt. ostflatgmax(ilmda)) then
         scale = exp((ostflatgmax(ilmda)-g)/rt)
         lmdafsum(ilmda) = lmdafsum(ilmda) * scale
         lmdafwt(ilmda) = lmdafwt(ilmda) * scale
         ostflatgmax(ilmda) = g
      else
         w = exp((g-ostflatgmax(ilmda))/rt)
      end if
c
c     add the weighted sample and update the mean force of the bin
c
      lmdafsum(ilmda) = lmdafsum(ilmda) + w*lmdafhist(ihist)
      lmdafwt(ilmda) = lmdafwt(ilmda) + w
      lmdafmean(ilmda) = lmdafsum(ilmda) / lmdafwt(ilmda)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine clearkernels -- zero the ost bias kernels  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "clearkernels" empties the ost g kernel, its derivative grids
c     and running maxima, and the mean force of every lambda bin with
c     the reference bias of the sample mean force, first allocating
c     that reference if it is missing or no longer matches the number
c     of lambda bins
c
c     note the mean force of the lambda bins holds the kernel or the
c     sample mean force, whichever OST-FLATTEN selects, and never both,
c     so a change of the selection must be followed by this routine
c     and a rebuild of the kernels from the history
c
c
      subroutine clearkernels
      use dlmda
      use ost
      implicit none
      integer i,j
      logical resize
c
c
c     allocate the reference bias of the sample mean force if needed
c
      resize = (.not. allocated(ostflatgmax))
      if (.not. resize)  resize = (size(ostflatgmax) .ne. nlmda)
      if (resize) then
         if (allocated(ostflatgmax))  deallocate (ostflatgmax)
         allocate (ostflatgmax(nlmda))
      end if
c
c     empty the kernels, the mean force and the reference bias
c
      do i = 1, nlmda
         ostflatgmax(i) = 0.0d0
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
c     bias history; an existing history is read back and continued
c     at an interval boundary, and every setting stored in its header
c     must match the current run, while otherwise a new file is started
c
c
      subroutine initostfile
      use bath
      use dlmda
      use files
      use iounit
      use ost
      implicit none
      integer ihis
      integer nlmda0
      integer lmdaintv0
      integer freeunit
      integer trimtext
      real*8 wflmda0
      real*8 oststdev0,kelvin0
      real*8 lmdakelvin0
      real*8 lmdatemp
      real*8 lmdamass0,lmdafric0
      real*8 lmdadt0
      logical exist
      character*11 rname
c
c
c     return if ost has not been initialized
c
      if (.not. use_ost)  return
c
c     samples taken while lambda moves follow the joint distribution
c     of lambda and the atoms, which gives the mean force only when
c     the lambda particle moves at the temperature of the atoms
c
      if (.not. use_lmdahold) then
         if (abs(lmdatemp()-kelvin) .gt. 1.0d-10*max(1.0d0,kelvin)) then
            write (iout,5)
    5       format (/,' INITOSTFILE  --  LAMBDA-SAMPLE MOVING',
     &                 ' Requires LAMBDA-TEMPERATURE',
     &              /,' at the Bath Temperature')
            call fatal
         end if
      end if
c
c     continue from the latest existing ost history file
c
      lmdasavefile = filename(1:leng)//'.ost'
      call version (lmdasavefile,'old')
      inquire (file=lmdasavefile,exist=exist)
      if (exist) then
         nlmda0 = nlmda
         lmdaintv0 = lmdaintv
         wflmda0 = wflmda
         oststdev0 = oststdev
         kelvin0 = kelvin
         lmdakelvin0 = lmdatemp ()
         lmdamass0 = lmdamass
         lmdafric0 = lmdafric
         lmdadt0 = lmdadt
         call rdost
c
c     every setting stored in the header must match the current run
c
         rname = 'INITOSTFILE'
         call chkbiasset (rname,'LAMBDA-NBIN',dble(nlmda0),
     &                    dble(nlmda),lmdasavefile)
         call chkbiasset (rname,'LAMBDA-INTERVAL',dble(lmdaintv0),
     &                    dble(lmdaintv),lmdasavefile)
         call chkbiasset (rname,'FLAMBDA-WIDTH',wflmda0,wflmda,
     &                    lmdasavefile)
         call chkbiasset (rname,'OST-STDDEV',oststdev0,oststdev,
     &                    lmdasavefile)
         call chkbiasset (rname,'Temperature',kelvin0,kelvin,
     &                    lmdasavefile)
         call chkbiasset (rname,'LAMBDA-MASS',lmdamass0,lmdamass,
     &                    lmdasavefile)
         call chkbiasset (rname,'LAMBDA-TEMPERATURE',lmdakelvin0,
     &                    lmdakelvin,lmdasavefile)
         call chkbiasset (rname,'LAMBDA-FRICTION',lmdafric0,lmdafric,
     &                    lmdasavefile)
         call chkbiasset (rname,'LAMBDA-DT',lmdadt0,lmdadt,
     &                    lmdasavefile)
c
c     resume at an interval boundary, since the samples of a partial
c     interval are not saved in the history file
c
         lmdastep = (lmdastep/lmdaintv) * lmdaintv
c
c     drop the rows an interrupted save appended after the last
c     header update, or widen a header written before the lambda
c     particle temperature was stored, by rewriting the file from
c     the gaussians just read; a clean file is left alone until the
c     first save updates its header
c
         if (histrewrite) then
            ihis = freeunit ()
            open (unit=ihis,file=lmdasavefile,status='replace',
     &            access='stream',form='formatted')
            call prtost (ihis)
            close (unit=ihis)
         end if
         write (iout,10)  lmdasavefile(1:trimtext(lmdasavefile))
   10    format (/,' OST  --  Bias History Continued From  ',a)
         return
      end if
c
c     otherwise open a new file holding the current history
c
      ihis = freeunit ()
      call version (lmdasavefile,'new')
      open (unit=ihis,file=lmdasavefile,status='new',
     &      access='stream',form='formatted')
      call prtost (ihis)
      close (unit=ihis)
      nlmdasave = nlmdahist
c
c     report the name of the file holding the bias history
c
      write (iout,20)  lmdasavefile(1:trimtext(lmdasavefile))
   20 format (/,' OST  --  Bias History Written To  ',a)
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
c     kernels; used to continue an ost simulation
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
      real*8 efreetot
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
c     read the gaussian history entries
c
      call rdosthist (ihis,ostfile)
c
c     rebuild the kernels from the saved gaussians
c
      if (fastkernel) then
         call buildkernels
      else
         call buildgkernel
         call buildfkernel
      end if
      lmdadeltag = efreetot()
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
c     the histogram lookup; the caller builds the bias kernels
c
c
      subroutine rdosthist (ihis,ostfile)
      use dlmda
      use inform
      use iounit
      use ost
      implicit none
      integer i,k,ihis
      integer ihist
      integer osthist0
      integer trimtext
      logical fixed
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
c     initialize ost arrays, where the lookup table is cleared when
c     it is rebuilt below
c
      call clearkernels
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
c     read saved gaussian histogram entries, by the fixed format
c     of "prtosthist" when a row has its exact layout, which is much
c     faster than a free-format read; each of the seven fields must
c     end in its last column, with a blank starting the next field,
c     since blanks inside a field are ignored by a fixed-format read
c
      do ihist = 1, nlmdahist
         read (ihis,10,err=90,end=90)  record
         fixed = (trimtext(record) .eq. 154)
         do i = 1, 6
            k = 12 * min(i,2) + 26 * max(i-2,0)
            if (record(k:k).eq.' ' .or. record(k+1:k+1).ne.' ')
     &         fixed = .false.
         end do
         if (fixed) then
            read (record,13,err=12,end=12)  lmdaihist(ihist),osthist0,
     &         lmdalhist(ihist),lmdafhist(ihist),osthhist(ihist),
     &         ostwlhist(ihist),ostwfhist(ihist)
   13       format (2i12,5d26.16)
            goto 14
         end if
   12    continue
         read (record,*,err=90,end=90)  lmdaihist(ihist),osthist0,
     &      lmdalhist(ihist),lmdafhist(ihist),osthhist(ihist),
     &      ostwlhist(ihist),ostwfhist(ihist)
   14    continue
         osthist(ihist) = osthist0
         maxwlhist = max(maxwlhist,ostwlhist(ihist))
         maxwfhist = max(maxwfhist,ostwfhist(ihist))
      end do
   10 format (a240)
c
c     note any rows past the saved count, which an interrupted save
c     leaves behind, beside an old header that "rdbiashead" noted
c
      read (ihis,10,err=15,end=15)  record
      histrewrite = .true.
   15 continue
      close (unit=ihis)
c
c     rebuild the lookup table from the saved gaussians, where the
c     caller builds the kernels
c
      call buildostindex
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
c     metadynamics bias history; an existing history is read back and
c     continued at an interval boundary, and every setting stored in
c     its header must match the current run, while otherwise a new
c     file is started
c
c
      subroutine initmetafile
      use bath
      use dlmda
      use files
      use iounit
      use ost
      implicit none
      integer ihis
      integer nlmda0
      integer lmdaintv0
      integer freeunit
      integer trimtext
      real*8 oststdev0,kelvin0
      real*8 lmdakelvin0
      real*8 lmdatemp
      real*8 lmdamass0,lmdafric0
      real*8 lmdadt0
      logical exist
      character*12 rname
c
c
c     return if metadynamics has not been initialized
c
      if (.not. use_meta)  return
c
c     a Monte Carlo trial accepts on the instantaneous energy change,
c     so a lambda particle hotter or colder than the atoms does not
c     sample the free energy at its own temperature, however slowly
c     it moves, and the bias would not converge to the free energy
c
      if (use_lmdamc) then
         if (abs(lmdatemp()-kelvin) .gt. 1.0d-10*max(1.0d0,kelvin)) then
            write (iout,5)
    5       format (/,' INITMETAFILE  --  Monte Carlo Lambda Moves',
     &                 ' in Metadynamics Require',
     &              /,' LAMBDA-TEMPERATURE at the Bath Temperature')
            call fatal
         end if
      end if
c
c     continue from the latest existing metadynamics history file
c
      metasavefile = filename(1:leng)//'.meta'
      call version (metasavefile,'old')
      inquire (file=metasavefile,exist=exist)
      if (exist) then
         nlmda0 = nlmda
         lmdaintv0 = lmdaintv
         oststdev0 = oststdev
         kelvin0 = kelvin
         lmdakelvin0 = lmdatemp ()
         lmdamass0 = lmdamass
         lmdafric0 = lmdafric
         lmdadt0 = lmdadt
         call rdmeta
c
c     every setting stored in the header must match the current run
c
         rname = 'INITMETAFILE'
         call chkbiasset (rname,'LAMBDA-NBIN',dble(nlmda0),
     &                    dble(nlmda),metasavefile)
         call chkbiasset (rname,'LAMBDA-INTERVAL',dble(lmdaintv0),
     &                    dble(lmdaintv),metasavefile)
         call chkbiasset (rname,'OST-STDDEV',oststdev0,oststdev,
     &                    metasavefile)
         call chkbiasset (rname,'Temperature',kelvin0,kelvin,
     &                    metasavefile)
         call chkbiasset (rname,'LAMBDA-MASS',lmdamass0,lmdamass,
     &                    metasavefile)
         call chkbiasset (rname,'LAMBDA-TEMPERATURE',lmdakelvin0,
     &                    lmdakelvin,metasavefile)
         call chkbiasset (rname,'LAMBDA-FRICTION',lmdafric0,lmdafric,
     &                    metasavefile)
         call chkbiasset (rname,'LAMBDA-DT',lmdadt0,lmdadt,
     &                    metasavefile)
c
c     resume at an interval boundary, since the samples of a partial
c     interval are not saved in the history file
c
         lmdastep = (lmdastep/lmdaintv) * lmdaintv
c
c     drop the rows an interrupted save appended after the last
c     header update, or widen a header written before the lambda
c     particle temperature was stored, by rewriting the file from
c     the gaussians just read; a clean file is left alone until the
c     first save updates its header
c
         if (histrewrite) then
            ihis = freeunit ()
            open (unit=ihis,file=metasavefile,status='replace',
     &            access='stream',form='formatted')
            call prtmeta (ihis)
            close (unit=ihis)
         end if
         write (iout,10)  metasavefile(1:trimtext(metasavefile))
   10    format (/,' METADYN  --  Bias History Continued From  ',a)
         return
      end if
c
c     otherwise open a new file holding the current history
c
      ihis = freeunit ()
      call version (metasavefile,'new')
      open (unit=ihis,file=metasavefile,status='new',
     &      access='stream',form='formatted')
      call prtmeta (ihis)
      close (unit=ihis)
      nmethistsave = nmetahist
c
c     report the name of the file holding the bias history
c
      write (iout,20)  metasavefile(1:trimtext(metasavefile))
   20 format (/,' METADYN  --  Bias History Written To  ',a)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine rdmeta  --  input metadynamics history  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "rdmeta" reads a saved metadynamics bias history from the
c     external .meta file, then rebuilds the bin lookup and the
c     bias grid; used to continue a metadynamics simulation
c
c
      subroutine rdmeta
      use iounit
      use ost
      implicit none
      integer ihis
      integer freeunit
      integer trimtext
      logical exist
      character*240 metafile
      character*240 title
c
c
c     return unless a metadynamics history file is present
c
      metafile = metasavefile
      inquire (file=metafile,exist=exist)
      if (.not. exist)  return
      ihis = freeunit ()
      open (unit=ihis,file=metafile,status='old')
      rewind (unit=ihis)
c
c     read the header, which must be that of a metadynamics history
c
      call rdbiashead (ihis,metafile,title)
      if (index(title,metatitle(1:trimtext(metatitle))) .eq. 0) then
         close (unit=ihis)
         write (iout,10)  metafile(1:trimtext(metafile))
   10    format (/,' RDMETA  --  File is not a Metadynamics History',
     &              ' File',
     &           /,'            File Name :  ',a)
         call fatal
      end if
c
c     read the gaussian history entries and rebuild the grid
c
      call rdmetahist (ihis,metafile)
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
c
c     note any rows past the saved count, which an interrupted save
c     leaves behind, beside an old header that "rdbiashead" noted
c
      read (ihis,10,err=15,end=15)  record
      histrewrite = .true.
   15 continue
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
