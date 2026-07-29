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
c     propagate the lambda particle are saved for a later eostdyn
c     or emetadyn call in the same step.
c
c
      subroutine eostbias
      use atoms
      use deriv
      use dlmda
      use energi
      use ost
      use virial
      implicit none
      integer i,j
      real*8 egbias,dgdl,dgdfl
      real*8 eostlmda,dfdl
      real*8 vbias,dvdl
c
c
c     refresh the unbiased dU/dlambda at the current configuration
c
      ostdedl = dedl
c
c     the metadynamics bias depends on lambda alone
c
      if (use_metadyn) then
         call emetabias (ostlambda,vbias,dvdl)
         esum = esum + vbias
         ostbvbias = vbias
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
      call efkernel (eostlmda,dfdl)
      esum = esum + egbias - eostlmda
      ostbvbias = egbias - eostlmda
      ostbdgdl = dgdl
      ostbdgdfl = dgdfl
      ostbdfdl = dfdl
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
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine eostdyn -- orthogonal space tempering algorithm  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "eostdyn" calculates free energies using the orthogonal space
c     tempering algorithm using biasing gaussians
c
c
      subroutine eostdyn
      use dlmda
      use ost
      implicit none
      integer k
      integer isamp,istep
      integer ilmda,iflmda
      integer lambdabin,flambdabin
      real*8 etotfkernel
c
c
c     increment iost step counter
c
      iost = iost + 1
c
c     build the effective lambda derivative from the bias terms
c     evaluated this step by eostbias, which also refreshed ostdedl
c
      ostdgdl = ostbdgdl + ostbdgdfl*d2edl2
      ostddgdl = ostbdfdl
      deffdl = ostdedl + ostdgdl - ostddgdl
c
c     save all values in the hist interval, but average only after
c     the requested equilibration fraction
c
      istep = mod(iost,iosthist)
      if (istep .eq. 0) then
         isamp = iosthist
      else
         isamp = istep
      end if
      ostllist(isamp) = ostlambda
      ostflist(isamp) = ostdedl
c
c     add a new histogram count every iosthist steps
c
      if (istep .eq. 0) then
         call avgstd (ostllist,ostnequil+1,ostnavg,
     &                ostlambdaavg,ostlambdastd)
         call avgstd (ostflist,ostnequil+1,ostnavg,
     &                ostdedlavg,ostdedlstd)
         ilmda = lambdabin(ostlambdaavg)
c
c     ensure histogram contains the unbiased dU/dlambda value
c
         maxwlhist = max(maxwlhist,wlhist)
         maxwfhist = max(maxwfhist,wfhist)
         call ensureflambda (ostdedlavg)
         iflmda = flambdabin(ostdedlavg)
c
c     ensure histogram array is sufficiently large
c
         nosthist = nosthist + 1
         if (nosthist .gt. sizeosthist)  call resizeosthist
         call ij_to_k(ilmda,iflmda,nlmda,k)
c
c     save histogram information
c
         osthist(nosthist) = k
         ostihist(nosthist) = iost
         ostlhist(nosthist) = ostlambdaavg
         ostfhist(nosthist) = ostdedlavg
         osthhist(nosthist) = hbias
         ostwlhist(nosthist) = wlhist
         ostwfhist(nosthist) = wfhist
         ostnext(nosthist) = osthead(ilmda,iflmda)
         osthead(ilmda,iflmda) = nosthist
         if (fastkernel) then
            call updatekernels
         else
            call updategkernel
            call buildfkernel
         end if
         eosttot = etotfkernel()
      end if
c
c     propagate the lambda particle for the next dynamics step
c
      call ostlangevin
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine emetadyn -- 1D lambda metadynamics method  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "emetadyn" applies a one-dimensional metadynamics bias along
c     the main lambda coordinate and deposits lambda gaussians
c
c
      subroutine emetadyn
      use ost
      implicit none
      integer istep,isamp
      real*8 metadeltag
c
c
c     increment adaptive-bias step counter
c
      iost = iost + 1
c
c     effective lambda derivative from the bias evaluated this step
c     by eostbias, which also refreshed ostdedl
c
      deffdl = ostdedl + ostbdgdl
c
c     save all lambda values in the hist interval
c
      istep = mod(iost,iosthist)
      if (istep .eq. 0) then
         isamp = iosthist
      else
         isamp = istep
      end if
      ostllist(isamp) = ostlambda
c
c     add a new metadynamics gaussian every iosthist steps
c
      if (istep .eq. 0) then
         call avgstd (ostllist,ostnequil+1,ostnavg,
     &                ostlambdaavg,ostlambdastd)
         nmetahist = nmetahist + 1
         if (nmetahist .gt. sizemetahist)  call resizemeta
         metalhist(nmetahist) = ostlambdaavg
         metahhist(nmetahist) = hbias
         metawhist(nmetahist) = wlmda
         metaihist(nmetahist) = iost
         eosttot = metadeltag()
      end if
c
c     propagate the lambda particle for the next dynamics step
c
      call ostlangevin
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
c     "emetabias" evaluates Vbias(lambda) and dVbias/dlambda for
c     the sum of one-dimensional normalized metadynamics gaussians
c
c
      subroutine emetabias (lambda,vbias,dvdl)
      use math
      use ost
      implicit none
      integer ihist
      real*8 lambda
      real*8 vbias,dvdl
      real*8 delta
      real*8 sig,sig2
      real*8 bias,pref
c
c
c     initialize metadynamics bias and derivative
c
      vbias = 0.0d0
      dvdl = 0.0d0
c
c     loop over saved metadynamics gaussians
c
      do ihist = 1, nmetahist
         sig = metawhist(ihist)
         if (sig .gt. 0.0d0) then
            sig2 = sig * sig
            delta = lambda - metalhist(ihist)
            pref = metahhist(ihist) / (sig*sqrt(2.0d0*pi))
            bias = pref * exp(-0.5d0*delta*delta/sig2)
            vbias = vbias + bias
            dvdl = dvdl - delta*bias/sig2
         end if
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
c
      function metadeltag ()
      implicit none
      real*8 metadeltag
      real*8 vbias0,vbias1
      real*8 dvdl
c
c
c     evaluate endpoint metadynamics bias values
c
      call emetabias (0.0d0,vbias0,dvdl)
      call emetabias (1.0d0,vbias1,dvdl)
      metadeltag = -vbias1 + vbias0
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
      integer i
      integer oldsize,newsize
      integer, allocatable :: metaihist0(:)
      real*8, allocatable :: metalhist0(:)
      real*8, allocatable :: metahhist0(:)
      real*8, allocatable :: metawhist0(:)
c
c
c     save old metadynamics history
c
      oldsize = sizemetahist
      newsize = 2 * oldsize
      allocate (metalhist0(oldsize))
      allocate (metahhist0(oldsize))
      allocate (metawhist0(oldsize))
      allocate (metaihist0(oldsize))
      do i = 1, oldsize
         metalhist0(i) = metalhist(i)
         metahhist0(i) = metahhist(i)
         metawhist0(i) = metawhist(i)
         metaihist0(i) = metaihist(i)
      end do
c
c     allocate resized metadynamics history
c
      deallocate (metalhist)
      deallocate (metahhist)
      deallocate (metawhist)
      deallocate (metaihist)
      sizemetahist = newsize
      allocate (metalhist(sizemetahist))
      allocate (metahhist(sizemetahist))
      allocate (metawhist(sizemetahist))
      allocate (metaihist(sizemetahist))
      do i = 1, oldsize
         metalhist(i) = metalhist0(i)
         metahhist(i) = metahhist0(i)
         metawhist(i) = metawhist0(i)
         metaihist(i) = metaihist0(i)
      end do
      do i = oldsize+1, newsize
         metalhist(i) = 0.0d0
         metahhist(i) = 0.0d0
         metawhist(i) = 0.0d0
         metaihist(i) = 0
      end do
      deallocate (metalhist0)
      deallocate (metahhist0)
      deallocate (metawhist0)
      deallocate (metaihist0)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ostlangevin -- propagate ost lambda particle  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ostlangevin" propagates the auxiliary lambda particle in
c     theta space, where lambda = sin(theta)**2
c
c
      subroutine ostlangevin
      use bath
      use math
      use ost
      use units
      implicit none
      real*8 c
      real*8 force
      real*8 gamma
      real*8 normal
      real*8 sigma
      real*8 sinth
      real*8 ktm
      external normal
c
c
c     return if lambda dynamics parameters are invalid
c
      if (ostdt .le. 0.0d0)  return
      if (ostmass .le. 0.0d0)  return
c
c     force on theta from dU/dlambda and lambda = sin(theta)**2
c
      force = -deffdl * sin(2.0d0*osttheta)
c
c     propagate theta velocity with Langevin friction and noise
c
      gamma = max(0.0d0,ostfriction)
      if (gamma .gt. 0.0d0) then
         c = exp(-gamma*ostdt)
         ktm = boltzmann * kelvin / ostmass
         sigma = sqrt(ktm*(1.0d0-c*c))
         ostvtheta = c*ostvtheta
     &                 + (1.0d0-c)*force/(gamma*ostmass)
     &                 + sigma*normal()
      else
         ostvtheta = ostvtheta + ostdt*force/ostmass
      end if
c
c     update theta and wrap it into the principal periodic interval
c
      osttheta = osttheta + ostdt*ostvtheta
      do while (osttheta .gt. pi)
         osttheta = osttheta - 2.0d0*pi
      end do
      do while (osttheta .le. -pi)
         osttheta = osttheta + 2.0d0*pi
      end do
c
c     map theta back to lambda and refresh sublambda values
c
      sinth = sin(osttheta)
      ostlambda = sinth * sinth
      call mapsublmda (ostlambda)
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
      nfcut = int(oststdev*maxwfhist/wflmda)
      if (dble(nfcut)*wflmda .lt. oststdev*maxwfhist+wflmda2)
     &   nfcut = nfcut + 1
      nfcut = nfcut + nbuffer
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
      use math
      use ost
      implicit none
      integer ihist
      integer ilmda,iflmda
      integer klmda,kflmda
      integer lcenter,lcount
      integer flcenter
      integer nlcut,nfcut
      integer img,nimg
      integer lambdabin,flambdabin
      real*8 egbias,dgdl,dgdfl
      real*8 sigl,sigf
      real*8 sigl2,sigf2
      real*8 sigl2inv,sigf2inv
      real*8 pref
      real*8 sourcel,sourcefl
      real*8 ldelta,ldelta2
      real*8 fldelta,fldelta2
      real*8 bias,expl,expfl
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
      ilmda = lambdabin(ostlambda)
      iflmda = nint(ostdedl / wflmda) + fli0
      if (iflmda .lt. 1 .or. iflmda .gt. nflmda)  return
c
c     use max gaussian widths for conservative bin cutoffs
c
      nlcut = int(oststdev*maxwlhist/wlmda)
      if (dble(nlcut)*wlmda .lt. oststdev*maxwlhist)
     &   nlcut = nlcut + 1
      nlcut = nlcut + 1
      nfcut = int(oststdev*maxwfhist/wflmda)
      if (dble(nfcut)*wflmda .lt. oststdev*maxwfhist)
     &   nfcut = nfcut + 1
      nfcut = nfcut + 1
c
c     loop over nearby lambda bins and their mirror images
c
      do klmda = -nlcut, nlcut
         lcenter = ilmda + klmda
         lcount = lcenter
         if (lcount .lt. 1) then
            lcount = 2 - lcount
         else if (lcount .gt. nlmda) then
            lcount = 2*nlmda - lcount
         end if
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
                     sourcefl = ostfhist(ihist)
                     fldelta = ostdedl - sourcefl
                     if (abs(fldelta) .le. oststdev*sigf) then
                        sigl2 = sigl * sigl
                        sigf2 = sigf * sigf
                        sigl2inv = 1.0d0 / sigl2
                        sigf2inv = 1.0d0 / sigf2
                        pref = osthhist(ihist)
     &                            / (2.0d0*pi*sigl*sigf)
                        fldelta2 = fldelta * fldelta
                        expfl = exp(-0.5d0*fldelta2*sigf2inv)
                        nimg = 1
                        if (lcenter .eq. 1 .or.
     &                      lcenter .eq. nlmda)  nimg = 2
                        do img = 1, nimg
                           if (lcenter .lt. 1) then
                              sourcel = -ostlhist(ihist)
                           else if (lcenter .gt. nlmda) then
                              sourcel = 2.0d0 - ostlhist(ihist)
                           else if (img .eq. 2 .and.
     &                              lcenter .eq. 1) then
                              sourcel = -ostlhist(ihist)
                           else if (img .eq. 2 .and.
     &                              lcenter .eq. nlmda) then
                              sourcel = 2.0d0 - ostlhist(ihist)
                           else
                              sourcel = ostlhist(ihist)
                           end if
                           ldelta = ostlambda - sourcel
                           if (abs(ldelta) .le. oststdev*sigl) then
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
      if (ostlambda .lt. 0.0d0 .or. ostlambda .gt. 1.0d0)
     &   return
      if (ostdedl .lt. flstart .or. ostdedl .gt. flend)  return
c
c     locate the lower-left grid point of the interpolation cell
c
      if (ostlambda .ge. 1.0d0) then
         il0 = nlmda - 1
      else
         il0 = int(ostlambda/wlmda) + 1
         il0 = max(1,min(il0,nlmda-1))
      end if
      if (ostdedl .ge. flend) then
         if0 = nflmda - 1
      else
         if0 = int((ostdedl-flstart)/wflmda) + 1
         if0 = max(1,min(if0,nflmda-1))
      end if
      l0 = dble(il0-1) * wlmda
      f0 = dble(if0-fli0) * wflmda
      x = (ostlambda-l0) / wlmda
      y = (ostdedl-f0) / wflmda
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
c     ######################################################
c     ##                                                  ##
c     ##  function lambdabin -- get bin index for lambda  ##
c     ##                                                  ##
c     ######################################################
c
c
c     "lambdabin" computes the bin index for ostlambda
c
c
      function lambdabin (lambda)
      use ost
      implicit none
      integer lambdabin
      real*8 lambda
c
c
c     set lambdabin value
c
      lambdabin = nint(lambda / wlmda) + 1
      if (lambdabin .lt. 1)  lambdabin = 1
      if (lambdabin .gt. nlmda)  lambdabin = nlmda
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
      use ost
      implicit none
      integer i
      integer oldsize
      integer newsize
      integer, allocatable :: osthist0(:)
      integer, allocatable :: ostihist0(:)
      integer, allocatable :: ostnext0(:)
      real*8, allocatable :: ostlhist0(:)
      real*8, allocatable :: ostfhist0(:)
      real*8, allocatable :: osthhist0(:)
      real*8, allocatable :: ostwlhist0(:)
      real*8, allocatable :: ostwfhist0(:)
c
c
c     save the old histogram data and set the new size
c
      oldsize = sizeosthist
      newsize = 2 * oldsize
      allocate (osthist0(oldsize))
      allocate (ostihist0(oldsize))
      allocate (ostnext0(oldsize))
      allocate (ostlhist0(oldsize))
      allocate (ostfhist0(oldsize))
      allocate (osthhist0(oldsize))
      allocate (ostwlhist0(oldsize))
      allocate (ostwfhist0(oldsize))
      do i = 1, oldsize
         osthist0(i) = osthist(i)
         ostihist0(i) = ostihist(i)
         ostnext0(i) = ostnext(i)
         ostlhist0(i) = ostlhist(i)
         ostfhist0(i) = ostfhist(i)
         osthhist0(i) = osthhist(i)
         ostwlhist0(i) = ostwlhist(i)
         ostwfhist0(i) = ostwfhist(i)
      end do
c
c     allocate and initialize the resized histogram arrays
c
      deallocate (osthist)
      deallocate (ostihist)
      deallocate (ostnext)
      deallocate (ostlhist)
      deallocate (ostfhist)
      deallocate (osthhist)
      deallocate (ostwlhist)
      deallocate (ostwfhist)
      sizeosthist = newsize
      allocate (osthist(sizeosthist))
      allocate (ostihist(sizeosthist))
      allocate (ostnext(sizeosthist))
      allocate (ostlhist(sizeosthist))
      allocate (ostfhist(sizeosthist))
      allocate (osthhist(sizeosthist))
      allocate (ostwlhist(sizeosthist))
      allocate (ostwfhist(sizeosthist))
c
c     restore old histogram data into the resized arrays
c
      do i = 1, oldsize
         osthist(i) = osthist0(i)
         ostihist(i) = ostihist0(i)
         ostnext(i) = ostnext0(i)
         ostlhist(i) = ostlhist0(i)
         ostfhist(i) = ostfhist0(i)
         osthhist(i) = osthhist0(i)
         ostwlhist(i) = ostwlhist0(i)
         ostwfhist(i) = ostwfhist0(i)
      end do
      do i = oldsize+1, newsize
         osthist(i) = 0
         ostihist(i) = 0
         ostnext(i) = 0
         ostlhist(i) = 0.0d0
         ostfhist(i) = 0.0d0
         osthhist(i) = 0.0d0
         ostwlhist(i) = 0.0d0
         ostwfhist(i) = 0.0d0
      end do
      deallocate (osthist0)
      deallocate (ostihist0)
      deallocate (ostnext0)
      deallocate (ostlhist0)
      deallocate (ostfhist0)
      deallocate (osthhist0)
      deallocate (ostwlhist0)
      deallocate (ostwfhist0)
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
      use ost
      implicit none
      integer i,j
      integer ihist
      integer k
      integer ilmda,iflmda
      integer lambdabin,flambdabin
c
c
c     clear the bin heads and next links
c
      do i = 1, nlmda
         do j = 1, nflmda
            osthead(i,j) = 0
         end do
      end do
      do i = 1, sizeosthist
         ostnext(i) = 0
      end do
c
c     insert each saved histogram entry into its lambda/flambda bin
c
      do ihist = 1, nosthist
         ilmda = lambdabin(ostlhist(ihist))
         iflmda = flambdabin(ostfhist(ihist))
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
      use ost
      implicit none
      integer i,j
      integer ihist
c
c
c     zero out g kernel
c
      do i = 1, nlmda
         do j = 1, nflmda
            gkernel(i,j) = 0.0d0
         end do
      end do
c
c     loop over saved histogram sources
c
      do ihist = 1, nosthist
         call addgkernelhist (ihist)
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
      use ost
      implicit none
      integer i,j
      integer ihist
c
c
c     zero out kernels and free energy accumulators
c
      do i = 1, nlmda
         fkernel(i) = 0.0d0
         fsumkernel(i) = 0.0d0
         pfkernel(i) = 0.0d0
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
      do ihist = 1, nosthist
         call addkernelhist (ihist)
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
      use ost
      implicit none
c
c
c     add the most recently saved histogram source
c
      if (nosthist .gt. 0)  call addgkernelhist (nosthist)
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
      use ost
      implicit none
c
c
c     add the most recently saved histogram source
c
      if (nosthist .gt. 0)  call addkernelhist (nosthist)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine addgkernelhist -- add one histogram source  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "addgkernelhist" spreads one saved histogram source to the
c     ost g kernel using normalized gaussian biasing functions
c
c
      subroutine addgkernelhist (ihist)
      use ost
      use math
      implicit none
      integer ihist
      integer k
      integer lsrc,fsrc
      integer ilmda,iflmda
      integer img,llog
      integer ilmda1,ilmda2
      integer iflmda1,iflmda2
      integer nlcut,nfcut
      real*8 sigl,sigf
      real*8 sigl2,sigf2
      real*8 pref
      real*8 sourcel,sourcefl
      real*8 targetl,targetf
      real*8 ldelta,ldelta2
      real*8 fldelta,fldelta2
      real*8 e,expl,expfl
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
      sourcefl = ostfhist(ihist)
c
c     include target bins out to the requested gaussian cutoff
c
      nlcut = int(oststdev*sigl/wlmda)
      if (dble(nlcut)*wlmda .lt. oststdev*sigl+wlmda2)
     &   nlcut = nlcut + 1
      nfcut = int(oststdev*sigf/wflmda)
      if (dble(nfcut)*wflmda .lt. oststdev*sigf+wflmda2)
     &   nfcut = nfcut + 1
      iflmda1 = max(1,fsrc-nfcut)
      iflmda2 = min(nflmda,fsrc+nfcut)
c
c     use the physical lambda source plus its two mirror images
c
      do img = 1, 3
         if (img .eq. 1) then
            llog = lsrc
            sourcel = ostlhist(ihist)
         else if (img .eq. 2) then
            llog = 2 - lsrc
            sourcel = -ostlhist(ihist)
         else
            llog = 2*nlmda - lsrc
            sourcel = 2.0d0 - ostlhist(ihist)
         end if
c
c     only target bins within the gaussian cutoff receive a contribution
c
         ilmda1 = max(1,llog-nlcut)
         ilmda2 = min(nlmda,llog+nlcut)
         if (ilmda1 .le. ilmda2) then
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
                  gkernel(ilmda,iflmda) = gkernel(ilmda,iflmda) + e
   20             continue
               end do
   10          continue
            end do
         end if
      end do
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
c     "addkernelhist" spreads one saved histogram source to the
c     ost g kernel and updates the f kernel accumulators
c
c
      subroutine addkernelhist (ihist)
      use ost
      use math
      implicit none
      integer ihist
      integer k
      integer lsrc,fsrc
      integer ilmda,iflmda
      integer img,llog
      integer ilmda1,ilmda2
      integer iflmda1,iflmda2
      integer nlcut,nfcut
      real*8 sigl,sigf
      real*8 sigl2,sigf2
      real*8 pref
      real*8 sourcel,sourcefl
      real*8 targetl,targetf
      real*8 ldelta,ldelta2
      real*8 fldelta,fldelta2
      real*8 e,expl,expfl
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
      sourcefl = ostfhist(ihist)
c
c     include target bins out to the requested gaussian cutoff
c
      nlcut = int(oststdev*sigl/wlmda)
      if (dble(nlcut)*wlmda .lt. oststdev*sigl+wlmda2)
     &   nlcut = nlcut + 1
      nfcut = int(oststdev*sigf/wflmda)
      if (dble(nfcut)*wflmda .lt. oststdev*sigf+wflmda2)
     &   nfcut = nfcut + 1
      iflmda1 = max(1,fsrc-nfcut)
      iflmda2 = min(nflmda,fsrc+nfcut)
c
c     use the physical lambda source plus its two mirror images
c
      do img = 1, 3
         if (img .eq. 1) then
            llog = lsrc
            sourcel = ostlhist(ihist)
         else if (img .eq. 2) then
            llog = 2 - lsrc
            sourcel = -ostlhist(ihist)
         else
            llog = 2*nlmda - lsrc
            sourcel = 2.0d0 - ostlhist(ihist)
         end if
c
c     only target bins within the gaussian cutoff receive a contribution
c
         ilmda1 = max(1,llog-nlcut)
         ilmda2 = min(nlmda,llog+nlcut)
         if (ilmda1 .le. ilmda2) then
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
                  call addkernelpoint (ilmda,iflmda,e,ldelta,
     &                                 fldelta,sigl2,sigf2)
   20             continue
               end do
   10          continue
            end do
         end if
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
c     f kernel numerator and partition-function accumulator
c
c
      subroutine addkernelpoint (ilmda,iflmda,e,ldelta,
     &                           fldelta,sigl2,sigf2)
      use bath
      use ost
      use units
      implicit none
      integer ilmda,iflmda
      real*8 e
      real*8 ldelta,fldelta
      real*8 sigl2,sigf2
      real*8 flmda
      real*8 oldg,newg
      real*8 oldweight,newweight
      real*8 delweight
      real*8 dgdl,dgdfl,d2gdlfl
c
c
c     update the g kernel and adjust the f kernel accumulators
c
      oldg = gkernel(ilmda,iflmda)
      if (oldg .eq. 0.0d0) then
         oldweight = 0.0d0
      else
         oldweight = exp(oldg/(gasconst*kelvin))
      end if
      newg = oldg + e
      newweight = exp(newg/(gasconst*kelvin))
      delweight = newweight - oldweight
      flmda = dble(iflmda-fli0) * wflmda
      dgdl = -ldelta * e / sigl2
      dgdfl = -fldelta * e / sigf2
      d2gdlfl = ldelta * fldelta * e / (sigl2*sigf2)
      gfkernel(ilmda,iflmda) = gfkernel(ilmda,iflmda) + dgdfl
      gkernel(ilmda,iflmda) = newg
      glfkernel(ilmda,iflmda) = glfkernel(ilmda,iflmda) + d2gdlfl
      glkernel(ilmda,iflmda) = glkernel(ilmda,iflmda) + dgdl
      fsumkernel(ilmda) = fsumkernel(ilmda) + flmda*delweight
      pfkernel(ilmda) = pfkernel(ilmda) + delweight
      if (pfkernel(ilmda) .eq. 0.0d0) then
         fkernel(ilmda) = 0.0d0
      else
         fkernel(ilmda) = fsumkernel(ilmda) / pfkernel(ilmda)
      end if
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
      use ost
      use units
      implicit none
      integer ilmda,iflmda
      real*8 avgflambda
      real*8 partfunc
      real*8 flmda
      real*8 weight
c
c
c     loop over lambda bins
c
      do ilmda = 1, nlmda
         avgflambda = 0.0d0
         partfunc = 0.0d0
c
c     loop over flambda bins with nonzero g kernel support
c
         do iflmda = 1, nflmda
            if (gkernel(ilmda,iflmda) .ne. 0.0d0) then
               flmda = dble(iflmda-fli0) * wflmda
               weight = exp(gkernel(ilmda,iflmda)/(gasconst*kelvin))
               avgflambda = avgflambda + flmda*weight
               partfunc = partfunc + weight
            end if
         end do
c
c     set mean force for this lambda bin
c
         if (partfunc .eq. 0.0d0) then
            fkernel(ilmda) = 0.0d0
         else
            fkernel(ilmda) = avgflambda / partfunc
         end if
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function etotfkernel -- compute free energy from fkernel  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "etotfkernel" computes the total free energy change by
c     integrating the f kernel over lambda bins using the trapezoid
c     rule, with dDeltaG/dlambda = fkernel
c
c
      function etotfkernel ()
      use ost
      implicit none
      integer ilmda
      real*8 etotfkernel
c
c
c     initialize free energy
c
      etotfkernel = 0.0d0
c
c     integrate over lambda bins
c
      do ilmda = 1, nlmda-1
         etotfkernel = etotfkernel
     &                + 0.5d0*(fkernel(ilmda)+fkernel(ilmda+1))*wlmda
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine efkernel -- free energy and derivative from f  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "efkernel" computes DeltaG at the current ostlambda by
c     integrating the f kernel using linear interpolation, and also
c     returns dDeltaG/dlambda at the current ostlambda
c
c
      subroutine efkernel (eostlmda,dfdl)
      use ost
      implicit none
      integer ilmda0,ilmda1
      real*8 eostlmda,dfdl
      real*8 fl0,fl1
      real*8 lmda0,lmda1
      real*8 slope
      real*8 x
c
c
c     initialize free energy and derivative
c
      eostlmda = 0.0d0
      dfdl = 0.0d0
c
c     handle endpoint at lambda = 0
c
      if (ostlambda .le. 0.0d0) then
         dfdl = fkernel(1)
         return
      end if
c
c     integrate over lambda intervals
c
      do ilmda0 = 1, nlmda-1
         ilmda1 = ilmda0 + 1
         lmda0 = dble(ilmda0-1) * wlmda
         lmda1 = dble(ilmda1-1) * wlmda
         fl0 = fkernel(ilmda0)
         fl1 = fkernel(ilmda1)
         slope = (fl1-fl0) / wlmda
c
c     integrate only to ostlambda if it lies in this interval
c
         if (ostlambda .le. lmda1) then
            x = ostlambda - lmda0
            eostlmda = eostlmda + fl0*x + 0.5d0*slope*x*x
            dfdl = fl0 + slope*x
            return
         end if
c
c     otherwise integrate the full interval
c
         eostlmda = eostlmda
     &              + 0.5d0*(fl0+fl1)*wlmda
      end do
c
c     handle endpoint at lambda = 1
c
      dfdl = fkernel(nlmda)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine saveost -- output ost restart information  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "saveost" writes the information needed to restart an ost
c     simulation to the external .ost restart file; the fixed-size
c     header is updated in place and new history is appended
c
c
      subroutine saveost
      use dlmda
      use files
      use ost
      implicit none
      integer ihis
      integer freeunit
      logical exist
      character*240 ostfile
c
c
c     return if ost has not been initialized
c
      if (.not. use_ost)  return
      if (.not. allocated(osthist))  return
c
c     create a new restart file or append only new history entries;
c     append before updating the header so an interrupted write leaves
c     the previous history count as the last complete checkpoint
c
      ihis = freeunit ()
      ostfile = filename(1:leng)//'.ost'
      inquire (file=ostfile,exist=exist)
      if (.not. exist .or. nosthist.lt.nosthistsave) then
         open (unit=ihis,file=ostfile,status='replace',
     &         access='stream',form='formatted')
         call prtost (ihis)
         close (unit=ihis)
      else
         if (nosthist .gt. nosthistsave) then
            open (unit=ihis,file=ostfile,status='old',
     &            access='stream',form='formatted',position='append')
            call prtosthist (ihis,nosthistsave+1,nosthist)
            close (unit=ihis)
         end if
         open (unit=ihis,file=ostfile,status='old',
     &         access='stream',form='unformatted',position='rewind')
         call updosthead (ihis)
         close (unit=ihis)
      end if
      nosthistsave = nosthist
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine rdost -- input ost restart information  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "rdost" reads the information needed to restart an ost
c     simulation from the external .ost restart file, then rebuilds
c     the histogram lookup and bias kernels
c
c
      subroutine rdost
      use dlmda
      use files
      use inform
      use iounit
      use ost
      implicit none
      integer i,ihis
      integer ihist
      integer freeunit
      integer trimtext
      integer iost0,iosthist0
      integer nlmda0,nflmda0,fli00
      integer nosthist0,sizeosthist0
      integer osthist0
      real*8 wlmda0,wflmda0
      real*8 ostlambda0
      real*8 osttheta0,ostvtheta0
      real*8 ostmass0,ostfriction0,ostdt0
      real*8 eosttot0,oststdev0
      real*8 etotfkernel
      logical exist
      character*240 ostfile
      character*240 record
c
c
c     return unless an ost restart file is present
c
      ostfile = filename(1:leng)//'.ost'
      inquire (file=ostfile,exist=exist)
      if (.not. exist)  return
      ihis = freeunit ()
      open (unit=ihis,file=ostfile,status='old')
      rewind (unit=ihis)
c
c     read scalar ost restart state
c
      read (ihis,10,err=90,end=90)  record
      read (ihis,10,err=90,end=90)  record
      read (ihis,10,err=90,end=90)  record
      read (record,*,err=90,end=90)  iost0,iosthist0,
     &   nlmda0,nflmda0,fli00,nosthist0,sizeosthist0
      read (ihis,10,err=90,end=90)  record
      read (ihis,10,err=90,end=90)  record
      read (record,*,err=90,end=90)  wlmda0,wflmda0,oststdev0
      read (ihis,10,err=90,end=90)  record
      read (ihis,10,err=90,end=90)  record
      read (record,*,err=90,end=90)  ostlambda0,osttheta0,ostvtheta0,
     &   ostmass0,ostfriction0,ostdt0,eosttot0
   10 format (a240)
c
c     validate restart dimensions
c
      if (nlmda0 .lt. 2)  goto 90
      if (nflmda0 .lt. 1)  goto 90
      if (fli00 .lt. 1 .or. fli00 .gt. nflmda0)  goto 90
      if (iosthist0 .lt. 1)  goto 90
      if (nosthist0 .lt. 0)  goto 90
      if (sizeosthist0 .lt. nosthist0)  sizeosthist0 = nosthist0
      if (sizeosthist0 .lt. 1)  sizeosthist0 = 1
c
c     reallocate ost arrays to match the restart file
c
      if (allocated(osthhist))  deallocate (osthhist)
      if (allocated(osthist))  deallocate (osthist)
      if (allocated(ostihist))  deallocate (ostihist)
      if (allocated(osthead))  deallocate (osthead)
      if (allocated(ostnext))  deallocate (ostnext)
      if (allocated(ostllist))  deallocate (ostllist)
      if (allocated(ostflist))  deallocate (ostflist)
      if (allocated(ostlhist))  deallocate (ostlhist)
      if (allocated(ostfhist))  deallocate (ostfhist)
      if (allocated(ostwlhist))  deallocate (ostwlhist)
      if (allocated(ostwfhist))  deallocate (ostwfhist)
      if (allocated(fkernel))  deallocate (fkernel)
      if (allocated(fsumkernel))  deallocate (fsumkernel)
      if (allocated(gfkernel))  deallocate (gfkernel)
      if (allocated(gkernel))  deallocate (gkernel)
      if (allocated(glfkernel))  deallocate (glfkernel)
      if (allocated(glkernel))  deallocate (glkernel)
      if (allocated(pfkernel))  deallocate (pfkernel)
      allocate (osthhist(sizeosthist0))
      allocate (osthist(sizeosthist0))
      allocate (ostihist(sizeosthist0))
      allocate (osthead(nlmda0,nflmda0))
      allocate (ostnext(sizeosthist0))
      allocate (ostllist(iosthist0))
      allocate (ostflist(iosthist0))
      allocate (ostlhist(sizeosthist0))
      allocate (ostfhist(sizeosthist0))
      allocate (ostwlhist(sizeosthist0))
      allocate (ostwfhist(sizeosthist0))
      allocate (fkernel(nlmda0))
      allocate (fsumkernel(nlmda0))
      allocate (gfkernel(nlmda0,nflmda0))
      allocate (gkernel(nlmda0,nflmda0))
      allocate (glfkernel(nlmda0,nflmda0))
      allocate (glkernel(nlmda0,nflmda0))
      allocate (pfkernel(nlmda0))
c
c     set scalar ost state from the restart file
c
      iost = iost0
      iosthist = iosthist0
      ostnequil = int(osteqratio*dble(iosthist))
      ostnequil = max(0,min(ostnequil,iosthist-1))
      ostnavg = iosthist - ostnequil
      nlmda = nlmda0
      nflmda = nflmda0
      fli0 = fli00
      nosthist = nosthist0
      sizeosthist = sizeosthist0
      wlmda = wlmda0
      wflmda = wflmda0
      wlmda2 = 0.5d0 * wlmda
      wflmda2 = 0.5d0 * wflmda
      maxwlhist = wlhist
      maxwfhist = wfhist
      ostlambda = ostlambda0
      ostlambdaavg = 0.0d0
      ostlambdastd = 0.0d0
      ostdedlavg = 0.0d0
      ostdedlstd = 0.0d0
      osttheta = osttheta0
      ostvtheta = ostvtheta0
      ostmass = ostmass0
      ostfriction = ostfriction0
      ostdt = ostdt0
      eosttot = eosttot0
      oststdev = oststdev0
c
c     initialize ost arrays
c
      do i = 1, nlmda
         fkernel(i) = 0.0d0
         fsumkernel(i) = 0.0d0
         pfkernel(i) = 0.0d0
         do ihist = 1, nflmda
            gfkernel(i,ihist) = 0.0d0
            gkernel(i,ihist) = 0.0d0
            glfkernel(i,ihist) = 0.0d0
            glkernel(i,ihist) = 0.0d0
            osthead(i,ihist) = 0
         end do
      end do
      do i = 1, iosthist
         ostllist(i) = 0.0d0
         ostflist(i) = 0.0d0
      end do
      do i = 1, sizeosthist
         osthist(i) = 0
         ostihist(i) = 0
         ostnext(i) = 0
         ostlhist(i) = 0.0d0
         ostfhist(i) = 0.0d0
         osthhist(i) = 0.0d0
         ostwlhist(i) = 0.0d0
         ostwfhist(i) = 0.0d0
      end do
c     read saved gaussian histogram entries
c
      read (ihis,10,err=90,end=90)  record
      do ihist = 1, nosthist
         read (ihis,10,err=90,end=90)  record
         read (record,*,err=90,end=90)  ostihist(ihist),osthist0,
     &      ostlhist(ihist),ostfhist(ihist),osthhist(ihist),
     &      ostwlhist(ihist),ostwfhist(ihist)
         osthist(ihist) = osthist0
         maxwlhist = max(maxwlhist,ostwlhist(ihist))
         maxwfhist = max(maxwfhist,ostwfhist(ihist))
      end do
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
      eosttot = etotfkernel()
      ostrestart = .true.
      nosthistsave = nosthist
      if (debug) then
         write (iout,20)  ostfile(1:trimtext(ostfile))
   20    format (/,' Restarting OST Bias from :  ',a)
      end if
      return
c
c     malformed restart file
c
   90 continue
      close (unit=ihis)
      write (iout,30)  ostfile(1:trimtext(ostfile))
   30 format (/,' RDOST  --  Error while Reading OST Restart File',
     &        /,'            File Name :  ',a)
      use_ost = .false.
      ostrestart = .false.
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
c     "prtost" writes current ost restart state
c
c
      subroutine prtost (ihis)
      use ost
      implicit none
      integer ihis
c
c     write the header followed by all current history entries
c
      call prtosthead (ihis)
      call prtosthist (ihis,1,nosthist)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine prtosthead  --  output ost restart header  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "prtosthead" writes the fixed-size current ost restart header
c
c
      subroutine prtosthead (ihis)
      use ost
      implicit none
      integer ihis
c
c
c     write scalar ost restart state
c
      write (ihis,10)
      write (ihis,20)
      write (ihis,30)  iost,iosthist,nlmda,nflmda,
     &                 fli0,nosthist,sizeosthist
      write (ihis,40)
      write (ihis,50)  wlmda,wflmda,oststdev
      write (ihis,60)
      write (ihis,70)  ostlambda,osttheta,ostvtheta,
     &                 ostmass,ostfriction,ostdt,eosttot
      write (ihis,80)
   10 format (' Orthogonal Space Tempering Restart :')
   20 format (' Integer State :')
   30 format (7i12)
   40 format (' Grid State :')
   50 format (3d26.16)
   60 format (' Lambda State :')
   70 format (7d26.16)
   80 format (' Gaussian History :')
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine updosthead  --  update ost restart header  ##
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
      integer ieol
      integer leol
      character*240 record
      character*2 newline
c
c     preserve the file's existing line-ending convention
c
      read (ihis,pos=1)  record
      ieol = index(record,achar(10))
      newline = achar(10)//' '
      leol = 1
      if (ieol .gt. 1) then
         if (record(ieol-1:ieol-1) .eq. achar(13)) then
            newline = achar(13)//achar(10)
            leol = 2
         end if
      end if
c
c     format each header record internally and write its raw bytes
c
      write (record,10)
      write (ihis,pos=1)  record(1:len_trim(record)),newline(1:leol)
      write (record,20)
      write (ihis)  record(1:len_trim(record)),newline(1:leol)
      write (record,30)  iost,iosthist,nlmda,nflmda,
     &                   fli0,nosthist,sizeosthist
      write (ihis)  record(1:len_trim(record)),newline(1:leol)
      write (record,40)
      write (ihis)  record(1:len_trim(record)),newline(1:leol)
      write (record,50)  wlmda,wflmda,oststdev
      write (ihis)  record(1:len_trim(record)),newline(1:leol)
      write (record,60)
      write (ihis)  record(1:len_trim(record)),newline(1:leol)
      write (record,70)  ostlambda,osttheta,ostvtheta,
     &                   ostmass,ostfriction,ostdt,eosttot
      write (ihis)  record(1:len_trim(record)),newline(1:leol)
      write (record,80)
      write (ihis)  record(1:len_trim(record)),newline(1:leol)
   10 format (' Orthogonal Space Tempering Restart :')
   20 format (' Integer State :')
   30 format (7i12)
   40 format (' Grid State :')
   50 format (3d26.16)
   60 format (' Lambda State :')
   70 format (7d26.16)
   80 format (' Gaussian History :')
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
         write (ihis,10)  ostihist(ihist),osthist(ihist),
     &      ostlhist(ihist),ostfhist(ihist),osthhist(ihist),
     &      ostwlhist(ihist),ostwfhist(ihist)
      end do
   10 format (2i12,5d26.16)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine savemeta -- output metadynamics restart  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "savemeta" writes the information needed to restart a
c     metadynamics simulation to the external .meta restart file; the
c     fixed-size header is updated in place and new history is appended
c
c
      subroutine savemeta
      use dlmda
      use files
      use ost
      implicit none
      integer ihis
      integer freeunit
      logical exist
      character*240 metafile
c
c
c     return if metadynamics has not been initialized
c
      if (.not. use_meta)  return
      if (.not. allocated(metalhist))  return
c
c     create a new restart file or append only new history entries
c
      ihis = freeunit ()
      metafile = filename(1:leng)//'.meta'
      inquire (file=metafile,exist=exist)
      if (.not. exist .or. nmetahist.lt.nmethistsave) then
         open (unit=ihis,file=metafile,status='replace',
     &         access='stream',form='formatted')
         call prtmeta (ihis)
         close (unit=ihis)
      else
         if (nmetahist .gt. nmethistsave) then
            open (unit=ihis,file=metafile,status='old',
     &            access='stream',form='formatted',position='append')
            call prtmetahist (ihis,nmethistsave+1,nmetahist)
            close (unit=ihis)
         end if
         open (unit=ihis,file=metafile,status='old',
     &         access='stream',form='unformatted',position='rewind')
         call updmetahead (ihis)
         close (unit=ihis)
      end if
      nmethistsave = nmetahist
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine rdmeta -- input metadynamics restart  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "rdmeta" reads the information needed to restart a
c     metadynamics simulation from the external .meta restart file
c
c
      subroutine rdmeta
      use dlmda
      use files
      use inform
      use iounit
      use ost
      implicit none
      integer i,ihis
      integer ihist
      integer freeunit
      integer trimtext
      integer iost0,iosthist0
      integer nmetahist0,sizemetahist0
      real*8 wlmda0
      real*8 ostlambda0
      real*8 osttheta0,ostvtheta0
      real*8 ostmass0,ostfriction0,ostdt0
      real*8 eosttot0
      real*8 metadeltag
      logical exist
      character*240 metafile
      character*240 record
c
c
c     return unless a metadynamics restart file is present
c
      metafile = filename(1:leng)//'.meta'
      inquire (file=metafile,exist=exist)
      if (.not. exist)  return
      ihis = freeunit ()
      open (unit=ihis,file=metafile,status='old')
      rewind (unit=ihis)
c
c     read scalar metadynamics restart state
c
      read (ihis,10,err=90,end=90)  record
      read (ihis,10,err=90,end=90)  record
      read (ihis,10,err=90,end=90)  record
      read (record,*,err=90,end=90)  iost0,iosthist0,
     &   nmetahist0,sizemetahist0
      read (ihis,10,err=90,end=90)  record
      read (ihis,10,err=90,end=90)  record
      read (record,*,err=90,end=90)  wlmda0
      read (ihis,10,err=90,end=90)  record
      read (ihis,10,err=90,end=90)  record
      read (record,*,err=90,end=90)  ostlambda0,osttheta0,ostvtheta0,
     &   ostmass0,ostfriction0,ostdt0,eosttot0
   10 format (a240)
c
c     validate restart dimensions
c
      if (iosthist0 .lt. 1)  goto 90
      if (nmetahist0 .lt. 0)  goto 90
      if (sizemetahist0 .lt. nmetahist0)
     &   sizemetahist0 = nmetahist0
      if (sizemetahist0 .lt. 1)  sizemetahist0 = 1
c
c     reallocate metadynamics arrays to match the restart file
c
      if (allocated(metalhist))  deallocate (metalhist)
      if (allocated(metahhist))  deallocate (metahhist)
      if (allocated(metawhist))  deallocate (metawhist)
      if (allocated(metaihist))  deallocate (metaihist)
      if (allocated(ostllist))  deallocate (ostllist)
      allocate (metalhist(sizemetahist0))
      allocate (metahhist(sizemetahist0))
      allocate (metawhist(sizemetahist0))
      allocate (metaihist(sizemetahist0))
      allocate (ostllist(iosthist0))
c
c     set scalar metadynamics state from the restart file
c
      iost = iost0
      iosthist = iosthist0
      ostnequil = int(osteqratio*dble(iosthist))
      ostnequil = max(0,min(ostnequil,iosthist-1))
      ostnavg = iosthist - ostnequil
      nmetahist = nmetahist0
      sizemetahist = sizemetahist0
      wlmda = wlmda0
      wlmda2 = 0.5d0 * wlmda
      ostlambda = ostlambda0
      ostlambdaavg = 0.0d0
      osttheta = osttheta0
      ostvtheta = ostvtheta0
      ostmass = ostmass0
      ostfriction = ostfriction0
      ostdt = ostdt0
      eosttot = eosttot0
c
c     initialize metadynamics arrays
c
      do i = 1, sizemetahist
         metalhist(i) = 0.0d0
         metahhist(i) = 0.0d0
         metawhist(i) = 0.0d0
         metaihist(i) = 0
      end do
      do i = 1, iosthist
         ostllist(i) = 0.0d0
      end do
c
c     read saved metadynamics gaussian entries
c
      read (ihis,10,err=90,end=90)  record
      do ihist = 1, nmetahist
         read (ihis,10,err=90,end=90)  record
         read (record,*,err=90,end=90)  metaihist(ihist),
     &      metalhist(ihist),metahhist(ihist),metawhist(ihist)
      end do
      close (unit=ihis)
      eosttot = metadeltag()
      metarestart = .true.
      nmethistsave = nmetahist
      if (debug) then
         write (iout,20)  metafile(1:trimtext(metafile))
   20    format (/,' Restarting Metadynamics Bias from :  ',a)
      end if
      return
c
c     malformed restart file
c
   90 continue
      close (unit=ihis)
      write (iout,30)  metafile(1:trimtext(metafile))
   30 format (/,' RDMETA  --  Error while Reading Metadynamics',
     &           ' Restart File',/,'             File Name :  ',a)
      use_meta = .false.
      metarestart = .false.
      call fatal
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine prtmeta -- output metadynamics restart  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "prtmeta" writes current metadynamics restart state
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
c     "prtmetahead" writes the fixed-size metadynamics restart header
c
c
      subroutine prtmetahead (ihis)
      use ost
      implicit none
      integer ihis
c
c
c     write scalar metadynamics restart state
c
      write (ihis,10)
      write (ihis,20)
      write (ihis,30)  iost,iosthist,nmetahist,sizemetahist
      write (ihis,40)
      write (ihis,50)  wlmda
      write (ihis,60)
      write (ihis,70)  ostlambda,osttheta,ostvtheta,
     &                 ostmass,ostfriction,ostdt,eosttot
      write (ihis,80)
   10 format (' Metadynamics Restart :')
   20 format (' Integer State :')
   30 format (4i12)
   40 format (' Grid State :')
   50 format (1d26.16)
   60 format (' Lambda State :')
   70 format (7d26.16)
   80 format (' Gaussian History :')
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
      integer ieol
      integer leol
      character*240 record
      character*2 newline
c
c     preserve the file's existing line-ending convention
c
      read (ihis,pos=1)  record
      ieol = index(record,achar(10))
      newline = achar(10)//' '
      leol = 1
      if (ieol .gt. 1) then
         if (record(ieol-1:ieol-1) .eq. achar(13)) then
            newline = achar(13)//achar(10)
            leol = 2
         end if
      end if
c
c     format each header record internally and write its raw bytes
c
      write (record,10)
      write (ihis,pos=1)  record(1:len_trim(record)),newline(1:leol)
      write (record,20)
      write (ihis)  record(1:len_trim(record)),newline(1:leol)
      write (record,30)  iost,iosthist,nmetahist,sizemetahist
      write (ihis)  record(1:len_trim(record)),newline(1:leol)
      write (record,40)
      write (ihis)  record(1:len_trim(record)),newline(1:leol)
      write (record,50)  wlmda
      write (ihis)  record(1:len_trim(record)),newline(1:leol)
      write (record,60)
      write (ihis)  record(1:len_trim(record)),newline(1:leol)
      write (record,70)  ostlambda,osttheta,ostvtheta,
     &                   ostmass,ostfriction,ostdt,eosttot
      write (ihis)  record(1:len_trim(record)),newline(1:leol)
      write (record,80)
      write (ihis)  record(1:len_trim(record)),newline(1:leol)
   10 format (' Metadynamics Restart :')
   20 format (' Integer State :')
   30 format (4i12)
   40 format (' Grid State :')
   50 format (1d26.16)
   60 format (' Lambda State :')
   70 format (7d26.16)
   80 format (' Gaussian History :')
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
