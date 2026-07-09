c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
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
      use atoms
      use deriv
      use dlmda
      use energi
      use ost
      use virial
      implicit none
      integer i,j
      integer k
      integer isamp,istep
      integer ilmda,iflmda
      integer lambdabin,flambdabin
      real*8 egbias,dgdl,dgdfl
      real*8 eostlmda,dfdl
      real*8 etotfkernel
c
c
c     increment iost step counter
c
      iost = iost + 1
c
c     compute current continuous g kernel bias and derivatives
c
      call egkernel (egbias,dgdl,dgdfl)
      call efkernel (eostlmda,dfdl)
      esum = esum + egbias - eostlmda
      ostdedl = dedl
      deffdl = ostdedl + dgdl + dgdfl*d2edl2 - dfdl
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
         call ostavgstd
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
         ostlhist(nosthist) = ostlambdaavg
         ostfhist(nosthist) = ostdedlavg
         osthhist(nosthist) = hbias
         ostwlhist(nosthist) = wlhist
         ostwfhist(nosthist) = wfhist
         ostnext(nosthist) = osthead(ilmda,iflmda)
         osthead(ilmda,iflmda) = nosthist
         call updategkernel
         call buildfkernel
         eosttot = etotfkernel()
c
c     zero out lambda and dE/dlambda averages and standard deviations
c
         ostlambdaavg = 0.0d0
         ostdedlavg = 0.0d0
         ostlambdastd = 0.0d0
         ostdedlstd = 0.0d0
      end if
c
c     propagate the lambda particle for the next dynamics step
c
      call ostlangevin
      return 
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ostavgstd -- ost average and std deviation    ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ostavgstd" computes average and standard deviation values
c     from lambda and dE/dlambda samples saved between hist updates
c
c
      subroutine ostavgstd
      use ost
      implicit none
      integer i
      real*8 dlmda,dfdl
c
c
c     compute averages from collected samples
c
      ostlambdaavg = 0.0d0
      ostdedlavg = 0.0d0
      do i = ostnequil+1, iosthist
         ostlambdaavg = ostlambdaavg + ostllist(i)
         ostdedlavg = ostdedlavg + ostflist(i)
      end do
      ostlambdaavg = ostlambdaavg / dble(ostnavg)
      ostdedlavg = ostdedlavg / dble(ostnavg)
c
c     compute population standard deviations from collected samples
c
      ostlambdastd = 0.0d0
      ostdedlstd = 0.0d0
      do i = ostnequil+1, iosthist
         dlmda = ostllist(i) - ostlambdaavg
         dfdl = ostflist(i) - ostdedlavg
         ostlambdastd = ostlambdastd + dlmda*dlmda
         ostdedlstd = ostdedlstd + dfdl*dfdl
      end do
      ostlambdastd = sqrt(ostlambdastd/dble(ostnavg))
      ostdedlstd = sqrt(ostdedlstd/dble(ostnavg))
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine emetadyn -- 1D lambda metadynamics method    ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "emetadyn" applies a one-dimensional metadynamics bias along
c     the main lambda coordinate and deposits lambda gaussians
c
c
      subroutine emetadyn
      use dlmda
      use energi
      use ost
      implicit none
      integer istep,navg
      real*8 vbias,dvdl
      real*8 metadeltag
c
c
c     increment adaptive-bias step counter
c
      iost = iost + 1
c
c     evaluate current metadynamics bias and add dVbias/dlambda
c
      call emetabias (ostlambda,vbias,dvdl)
      esum = esum + vbias
      ostdedl = dedl
      deffdl = ostdedl + dvdl
c
c     update average lambda only in the second half of each interval
c
      istep = mod(iost,iosthist)
      navg = iosthist / 2
      if ((istep .eq. 0) .or. (istep .gt. navg)) then
         ostlambdaavg = ostlambdaavg + ostlambda/dble(navg)
      end if
c
c     add a new metadynamics gaussian every iosthist steps
c
      if (istep .eq. 0) then
         nmetahist = nmetahist + 1
         if (nmetahist .gt. sizemetahist)  call resizemeta
         metalhist(nmetahist) = ostlambdaavg
         metahhist(nmetahist) = hbias
         metawhist(nmetahist) = wlmda
         eosttot = metadeltag()
         ostlambdaavg = 0.0d0
      end if
c
c     propagate the lambda particle for the next dynamics step
c
      call ostlangevin
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine emetabias -- evaluate metadynamics bias      ##
c     ##                                                           ##
c     ###############################################################
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
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine resizemeta -- resize metadynamics history     ##
c     ##                                                           ##
c     ###############################################################
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
      do i = 1, oldsize
         metalhist0(i) = metalhist(i)
         metahhist0(i) = metahhist(i)
         metawhist0(i) = metawhist(i)
      end do
c
c     allocate resized metadynamics history
c
      deallocate (metalhist)
      deallocate (metahhist)
      deallocate (metawhist)
      sizemetahist = newsize
      allocate (metalhist(sizemetahist))
      allocate (metahhist(sizemetahist))
      allocate (metawhist(sizemetahist))
      do i = 1, oldsize
         metalhist(i) = metalhist0(i)
         metahhist(i) = metahhist0(i)
         metawhist(i) = metawhist0(i)
      end do
      do i = oldsize+1, newsize
         metalhist(i) = 0.0d0
         metahhist(i) = 0.0d0
         metawhist(i) = 0.0d0
      end do
      deallocate (metalhist0)
      deallocate (metahhist0)
      deallocate (metawhist0)
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
      call mapsublmda
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine lmdachain -- chain rule for main lambda deriv  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "lmdachain" applies chain rule to the main lambda derivative
c     to compute the global lambda derivative of energy, energy^2,
c     force, and virial
c
c
      subroutine lmdachain
      use atoms
      use dlmda
      use mutant
      use ost
      implicit none
      integer i,j
c
c
c     apply chain rule for derivative of energy wrt global lambda
c
      d2epdl2 = d2epdl2 * dpldlmda*dpldlmda + depdl * d2pldlmda2
      depdl = depdl * dpldlmda
      d2evdl2 = d2evdl2 * dvldlmda*dvldlmda + devdl * d2vldlmda2
      devdl = devdl * dvldlmda
      d2emdl2 = d2emdl2 * deldlmda*deldlmda + demdl * d2eldlmda2
      demdl = demdl * deldlmda
      do i = 1, n
         do j = 1, 3
            dfpdl(j,i) = dfpdl(j,i) * dpldlmda
            dfmdl(j,i) = dfmdl(j,i) * deldlmda
            dfvdl(j,i) = dfvdl(j,i) * dvldlmda
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            depvirdl(j,i) = depvirdl(j,i) * dpldlmda
            demvirdl(j,i) = demvirdl(j,i) * deldlmda
            devvirdl(j,i) = devvirdl(j,i) * dvldlmda
         end do
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ensureflambda -- resize flambda histogram     ##
c     ##                                                           ##
c     ###############################################################
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
      real*8, allocatable :: gkernel0(:,:)
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
      allocate (gkernel0(nlmda,oldnflmda))
      do i = 1, nlmda
         do jflmda = 1, oldnflmda
            gkernel0(i,jflmda) = gkernel(i,jflmda)
         end do
      end do
      deallocate (gkernel)
      allocate (gkernel(nlmda,nflmda))
      do i = 1, nlmda
         do jflmda = 1, nflmda
            gkernel(i,jflmda) = 0.0d0
         end do
      end do
      do i = 1, nlmda
         do jflmda = 1, oldnflmda
            gkernel(i,jflmda+offset) = gkernel0(i,jflmda)
         end do
      end do
      deallocate (gkernel0)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine mapsublmda -- map from lambda to sublambda  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "mapsublmda" maps from main lambda to sublambdas
c
c
      subroutine mapsublmda
      use dlmda
      use mutant
      use ost
      implicit none
      real*8 taper
      real*8 dtaper
      real*8 d2taper
      character*6 mode
c
c
c     map from lambda to sublambdas
c
      if (ostpmap .eq. 'EXP') then
         call sublmdaexp (ostlambda,ostepexp,plambda,
     &                    dpldlmda,d2pldlmda2)
      else if (ostpmap .eq. 'INV') then
         call sublmdainvpower (ostlambda,ostinvepn,ostinvepeps,plambda,
     &                         dpldlmda,d2pldlmda2)
      else
         mode = 'OSTPOL'
         call sublmdataper (mode,ostlambda,taper,dtaper,d2taper)
         plambda = 1.0d0 - taper
         dpldlmda = -dtaper
         d2pldlmda2 = -d2taper
      end if
      if (ostemap .eq. 'EXP') then
         call sublmdaexp (ostlambda,ostemexp,elambda,
     &                    deldlmda,d2eldlmda2)
      else if (ostemap .eq. 'INV') then
         call sublmdainvpower (ostlambda,ostinvemn,ostinvemeps,elambda,
     &                         deldlmda,d2eldlmda2)
      else
         mode = 'OSTELE'
         call sublmdataper (mode,ostlambda,taper,dtaper,d2taper)
         elambda = 1.0d0 - taper
         deldlmda = -dtaper
         d2eldlmda2 = -d2taper
      end if
      if (ostvmap .eq. 'EXP') then
         call sublmdaexp (ostlambda,ostevexp,vlambda,
     &                    dvldlmda,d2vldlmda2)
      else if (ostvmap .eq. 'INV') then
         call sublmdainvpower (ostlambda,ostinvevn,ostinveveps,vlambda,
     &                         dvldlmda,d2vldlmda2)
      else
         mode = 'OSTVDW'
         call sublmdataper (mode,ostlambda,taper,dtaper,d2taper)
         vlambda = 1.0d0 - taper
         dvldlmda = -dtaper
         d2vldlmda2 = -d2taper
      end if
c
c     set flags to compute polarization lambda derivative
c
      use_pol4i = (plambda .le. ostplmda1)
      use_pol4f = (plambda .ge. ostplmda0)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine sublmdaexp -- exponential sublambda mapping   ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "sublmdaexp" maps from main lambda to a sublambda using a
c     power law and returns the first and second lambda derivatives
c
c
      subroutine sublmdaexp (x,exponent,lmda,dlmda,d2lmda)
      implicit none
      integer exponent
      real*8 x
      real*8 lmda
      real*8 dlmda
      real*8 d2lmda
      real*8 expnt
c
c
c     initialize values and handle endpoints explicitly
c
      expnt = dble(exponent)
      if (x .le. 0.0d0) then
         lmda = 0.0d0
         if (exponent .eq. 1) then
            dlmda = 1.0d0
            d2lmda = 0.0d0
         else if (exponent .eq. 2) then
            dlmda = 0.0d0
            d2lmda = 2.0d0
         else
            dlmda = 0.0d0
            d2lmda = 0.0d0
         end if
         return
      else if (x .ge. 1.0d0) then
         lmda = 1.0d0
         dlmda = expnt
         d2lmda = expnt * (expnt-1.0d0)
         return
      end if
c
c     compute lambda^exponent and its derivatives
c
      lmda = x**exponent
      dlmda = expnt * x**(exponent-1)
      if (exponent .eq. 1) then
         d2lmda = 0.0d0
      else
         d2lmda = expnt * (expnt-1.0d0) * x**(exponent-2)
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine sublmdainvpower -- inverse-power mapping      ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "sublmdainvpower" maps from main lambda to sublambda using a
c     shifted inverse-power law and returns first and second lambda
c     derivatives
c
c
      subroutine sublmdainvpower (x,n,eps,lmda,dlmda,d2lmda)
      implicit none
      integer n
      real*8 x
      real*8 eps
      real*8 lmda
      real*8 dlmda
      real*8 d2lmda
      real*8 xval
      real*8 shift
      real*8 power
      real*8 root0
      real*8 denom
      real*8 base
c
c
c     set a bounded coordinate and handle the identity map directly
c
      xval = x
      if (xval .lt. 0.0d0)  xval = 0.0d0
      if (xval .gt. 1.0d0)  xval = 1.0d0
      if (n .le. 1) then
         lmda = xval
         dlmda = 1.0d0
         d2lmda = 0.0d0
         return
      end if
c
c     compute normalized shifted inverse-power map
c
      shift = eps
      if (shift .le. 0.0d0)  shift = 0.1d0
      power = 1.0d0 / dble(n)
      root0 = shift**power
      denom = (1.0d0+shift)**power - root0
      base = xval + shift
      lmda = (base**power-root0) / denom
      dlmda = power * base**(power-1.0d0) / denom
      d2lmda = power * (power-1.0d0)
     &           * base**(power-2.0d0) / denom
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine sublmdataper -- tapers the sublambda  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "sublmdataper" tapers the mapping from main lambda to
c     sublambda at the endpoints
c
c
      subroutine sublmdataper (mode,x,taper,dtaper,d2taper)
      use shunt
      implicit none
      real*8 taper
      real*8 dtaper
      real*8 d2taper
      real*8 x,x2,x3
      real*8 x4,x5
      character*6 mode
c
c
c     get taper coefficients from existing Tinker switch routine
c
      call switch (mode)
c
c     return if outside switching window
c
      if (x .le. cut) then
         taper = 1.0d0
         dtaper = 0.0d0
         d2taper = 0.0d0
         return
      else if (x .ge. off) then
         taper = 0.0d0
         dtaper = 0.0d0
         d2taper = 0.0d0
         return
      end if
c
c     compute the quintic taper and derivative
c
      x2 = x*x
      x3 = x2*x
      x4 = x2*x2
      x5 = x2*x3
      taper = c5*x5 + c4*x4 + c3*x3 + c2*x2 + c1*x + c0
      dtaper = 5.0d0*c5*x4 + 4.0d0*c4*x3 + 3.0d0*c3*x2
     &            + 2.0d0*c2*x + c1
      d2taper = 20.0d0*c5*x3 + 12.0d0*c4*x2 + 6.0d0*c3*x
     &             + 2.0d0*c2
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
      allocate (ostnext0(oldsize))
      allocate (ostlhist0(oldsize))
      allocate (ostfhist0(oldsize))
      allocate (osthhist0(oldsize))
      allocate (ostwlhist0(oldsize))
      allocate (ostwfhist0(oldsize))
      do i = 1, oldsize
         osthist0(i) = osthist(i)
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
      deallocate (ostnext)
      deallocate (ostlhist)
      deallocate (ostfhist)
      deallocate (osthhist)
      deallocate (ostwlhist)
      deallocate (ostwfhist)
      sizeosthist = newsize
      allocate (osthist(sizeosthist))
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
         ostnext(i) = ostnext0(i)
         ostlhist(i) = ostlhist0(i)
         ostfhist(i) = ostfhist0(i)
         osthhist(i) = osthhist0(i)
         ostwlhist(i) = ostwlhist0(i)
         ostwfhist(i) = ostwfhist0(i)
      end do
      do i = oldsize+1, newsize
         osthist(i) = 0
         ostnext(i) = 0
         ostlhist(i) = 0.0d0
         ostfhist(i) = 0.0d0
         osthhist(i) = 0.0d0
         ostwlhist(i) = 0.0d0
         ostwfhist(i) = 0.0d0
      end do
      deallocate (osthist0)
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
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine buildostindex -- build ost histogram index    ##
c     ##                                                           ##
c     ###############################################################
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
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine addgkernelhist -- add one histogram source    ##
c     ##                                                           ##
c     ###############################################################
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
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine saveost -- output ost restart information     ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "saveost" writes the information needed to restart an ost
c     simulation to the external .ost restart file
c
c
      subroutine saveost
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
c     update an existing ost restart file or open a new one
c
      ihis = freeunit ()
      ostfile = filename(1:leng)//'.ost'
      inquire (file=ostfile,exist=exist)
      if (exist) then
         open (unit=ihis,file=ostfile,status='old')
         rewind (unit=ihis)
      else
         open (unit=ihis,file=ostfile,status='new')
      end if
      call prtost (ihis)
      close (unit=ihis)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine rdost -- input ost restart information        ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "rdost" reads the information needed to restart an ost
c     simulation from the external .ost restart file, then rebuilds
c     the histogram lookup and bias kernels
c
c
      subroutine rdost
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
      integer ihist0,osthist0
      real*8 wlmda0,wflmda0,wlmda20,wflmda20
      real*8 ostlambda0,ostlambdaavg0,ostdedlavg0
      real*8 osttheta0,ostvtheta0
      real*8 ostmass0,ostfriction0,ostdt0
      real*8 hbias0,eosttot0,oststdev0
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
      read (record,*,err=90,end=90)  wlmda0,wflmda0,
     &   wlmda20,wflmda20
      read (ihis,10,err=90,end=90)  record
      read (ihis,10,err=90,end=90)  record
      read (record,*,err=90,end=90)  ostlambda0,ostlambdaavg0,
     &   ostdedlavg0,osttheta0,ostvtheta0,ostmass0,ostfriction0,ostdt0
      read (ihis,10,err=90,end=90)  record
      read (ihis,10,err=90,end=90)  record
      read (record,*,err=90,end=90)  hbias0,eosttot0,oststdev0
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
      if (allocated(osthead))  deallocate (osthead)
      if (allocated(ostnext))  deallocate (ostnext)
      if (allocated(ostllist))  deallocate (ostllist)
      if (allocated(ostflist))  deallocate (ostflist)
      if (allocated(ostlhist))  deallocate (ostlhist)
      if (allocated(ostfhist))  deallocate (ostfhist)
      if (allocated(ostwlhist))  deallocate (ostwlhist)
      if (allocated(ostwfhist))  deallocate (ostwfhist)
      if (allocated(fkernel))  deallocate (fkernel)
      if (allocated(gkernel))  deallocate (gkernel)
      allocate (osthhist(sizeosthist0))
      allocate (osthist(sizeosthist0))
      allocate (osthead(nlmda0,nflmda0))
      allocate (ostnext(sizeosthist0))
      allocate (ostllist(iosthist0))
      allocate (ostflist(iosthist0))
      allocate (ostlhist(sizeosthist0))
      allocate (ostfhist(sizeosthist0))
      allocate (ostwlhist(sizeosthist0))
      allocate (ostwfhist(sizeosthist0))
      allocate (fkernel(nlmda0))
      allocate (gkernel(nlmda0,nflmda0))
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
      wlmda2 = wlmda20
      wflmda2 = wflmda20
      maxwlhist = wlhist
      maxwfhist = wfhist
      ostlambda = ostlambda0
      ostlambdaavg = ostlambdaavg0
      ostlambdastd = 0.0d0
      ostdedlavg = ostdedlavg0
      ostdedlstd = 0.0d0
      osttheta = osttheta0
      ostvtheta = ostvtheta0
      ostmass = ostmass0
      ostfriction = ostfriction0
      ostdt = ostdt0
      hbias = hbias0
      eosttot = eosttot0
      oststdev = oststdev0
c
c     initialize ost arrays
c
      do i = 1, nlmda
         fkernel(i) = 0.0d0
         do ihist = 1, nflmda
            gkernel(i,ihist) = 0.0d0
            osthead(i,ihist) = 0
         end do
      end do
      do i = 1, iosthist
         ostllist(i) = 0.0d0
         ostflist(i) = 0.0d0
      end do
      do i = 1, sizeosthist
         osthist(i) = 0
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
         read (record,*,err=90,end=90)  ihist0,osthist0,
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
      call buildgkernel
      call buildfkernel
      eosttot = etotfkernel()
      ostrestart = .true.
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
      integer ihist
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
      write (ihis,50)  wlmda,wflmda,wlmda2,wflmda2
      write (ihis,60)
      write (ihis,70)  ostlambda,ostlambdaavg,ostdedlavg,
     &                 osttheta,ostvtheta,ostmass,ostfriction,ostdt
      write (ihis,80)
      write (ihis,90)  hbias,eosttot,oststdev
      write (ihis,100)
c
c     write saved gaussian histogram entries
c
      do ihist = 1, nosthist
         write (ihis,110)  ihist,osthist(ihist),ostlhist(ihist),
     &      ostfhist(ihist),osthhist(ihist),ostwlhist(ihist),
     &      ostwfhist(ihist)
      end do
   10 format (' Orthogonal Space Tempering Restart :')
   20 format (' Integer State :')
   30 format (7i12)
   40 format (' Grid State :')
   50 format (4d26.16)
   60 format (' Lambda State :')
   70 format (8d26.16)
   80 format (' Bias State :')
   90 format (3d26.16)
  100 format (' Gaussian History :')
  110 format (2i12,5d26.16)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine savemeta -- output metadynamics restart       ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "savemeta" writes the information needed to restart a
c     metadynamics simulation to the external .meta restart file
c
c
      subroutine savemeta
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
c     update an existing metadynamics restart file or open a new one
c
      ihis = freeunit ()
      metafile = filename(1:leng)//'.meta'
      inquire (file=metafile,exist=exist)
      if (exist) then
         open (unit=ihis,file=metafile,status='old')
         rewind (unit=ihis)
      else
         open (unit=ihis,file=metafile,status='new')
      end if
      call prtmeta (ihis)
      close (unit=ihis)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine rdmeta -- input metadynamics restart          ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "rdmeta" reads the information needed to restart a
c     metadynamics simulation from the external .meta restart file
c
c
      subroutine rdmeta
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
      integer ihist0
      real*8 wlmda0,wlmda20
      real*8 ostlambda0,ostlambdaavg0
      real*8 osttheta0,ostvtheta0
      real*8 ostmass0,ostfriction0,ostdt0
      real*8 hbias0,eosttot0
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
      read (record,*,err=90,end=90)  wlmda0,wlmda20
      read (ihis,10,err=90,end=90)  record
      read (ihis,10,err=90,end=90)  record
      read (record,*,err=90,end=90)  ostlambda0,ostlambdaavg0,
     &   osttheta0,ostvtheta0,ostmass0,ostfriction0,ostdt0
      read (ihis,10,err=90,end=90)  record
      read (ihis,10,err=90,end=90)  record
      read (record,*,err=90,end=90)  hbias0,eosttot0
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
      allocate (metalhist(sizemetahist0))
      allocate (metahhist(sizemetahist0))
      allocate (metawhist(sizemetahist0))
c
c     set scalar metadynamics state from the restart file
c
      iost = iost0
      iosthist = iosthist0
      nmetahist = nmetahist0
      sizemetahist = sizemetahist0
      wlmda = wlmda0
      wlmda2 = wlmda20
      ostlambda = ostlambda0
      ostlambdaavg = ostlambdaavg0
      osttheta = osttheta0
      ostvtheta = ostvtheta0
      ostmass = ostmass0
      ostfriction = ostfriction0
      ostdt = ostdt0
      hbias = hbias0
      eosttot = eosttot0
c
c     initialize metadynamics arrays
c
      do i = 1, sizemetahist
         metalhist(i) = 0.0d0
         metahhist(i) = 0.0d0
         metawhist(i) = 0.0d0
      end do
c
c     read saved metadynamics gaussian entries
c
      read (ihis,10,err=90,end=90)  record
      do ihist = 1, nmetahist
         read (ihis,10,err=90,end=90)  record
         read (record,*,err=90,end=90)  ihist0,
     &      metalhist(ihist),metahhist(ihist),metawhist(ihist)
      end do
      close (unit=ihis)
      eosttot = metadeltag()
      metarestart = .true.
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
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine prtmeta -- output metadynamics restart        ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prtmeta" writes current metadynamics restart state
c
c
      subroutine prtmeta (ihis)
      use ost
      implicit none
      integer ihist
      integer ihis
c
c
c     write scalar metadynamics restart state
c
      write (ihis,10)
      write (ihis,20)
      write (ihis,30)  iost,iosthist,nmetahist,sizemetahist
      write (ihis,40)
      write (ihis,50)  wlmda,wlmda2
      write (ihis,60)
      write (ihis,70)  ostlambda,ostlambdaavg,osttheta,
     &                 ostvtheta,ostmass,ostfriction,ostdt
      write (ihis,80)
      write (ihis,90)  hbias,eosttot
      write (ihis,100)
c
c     write saved metadynamics gaussian entries
c
      do ihist = 1, nmetahist
         write (ihis,110)  ihist,metalhist(ihist),
     &      metahhist(ihist),metawhist(ihist)
      end do
   10 format (' Metadynamics Restart :')
   20 format (' Integer State :')
   30 format (4i12)
   40 format (' Grid State :')
   50 format (2d26.16)
   60 format (' Lambda State :')
   70 format (7d26.16)
   80 format (' Bias State :')
   90 format (2d26.16)
  100 format (' Gaussian History :')
  110 format (i12,3d26.16)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ij_to_k -- convert 2D index to 1D index       ##
c     ##                                                           ##
c     ###############################################################
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
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine k_to_ij -- convert 1D index to 2D index       ##
c     ##                                                           ##
c     ###############################################################
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
