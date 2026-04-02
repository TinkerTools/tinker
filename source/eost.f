c
c
c     ##################################################################
c     ##     COPYRIGHT (C) 2026 by  Moses Chung, Lianqing Zheng,      ##
c     ##        Jayvee Abella, Michael Schnieders, Pengyu Ren,        ##
c     ##                     Wei Yang, Jay Ponder                     ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine eost -- orthogonal space tempering algorithm  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "eost" calculates free energies using the orthogonal space
c     tempering algorithm using biasing gaussians
c
c
      subroutine eost
      use dlmda
      use mutant
      use ost
      implicit none
c
c
c     xxx
c
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
      depdl2 = depdl2 * dpldlmda*dpldlmda
     &           + depdl * d2pldlmda2
      depdl = depdl * dpldlmda
      devdl2 = devdl2 * dvldlmda*dvldlmda
     &           + devdl * d2vldlmda2
      devdl = devdl * dvldlmda
      demdl2 = demdl2 * deldlmda*deldlmda
     &           + demdl * d2eldlmda2
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
      real*8 x
      real*8 taper
      real*8 dtaper
      real*8 d2taper
      character*6 mode
c
c
c     map from lambda to sublambdas
c
      mode = 'OSTPOL'
      call sublmdataper (mode,ostlambda,taper,dtaper,d2taper)
      plambda = 1.0d0 - taper
      dpldlmda = -dtaper
      d2pldlmda2 = -d2taper
      mode = 'OSTELE'
      call sublmdataper (mode,ostlambda,taper,dtaper,d2taper)
      elambda = 1.0d0 - taper
      deldlmda = -dtaper
      d2eldlmda2 = -d2taper
      mode = 'OSTVDW'
      call sublmdataper (mode,ostlambda,taper,dtaper,d2taper)
      vlambda = 1.0d0 - taper
      dvldlmda = -dtaper
      d2vldlmda2 = -d2taper
c
c     set flags to compute polarization lambda derivative
c
      use_pol4i = (plambda .le. ostplmda1)
      use_pol4f = (plambda .ge. ostplmda0)
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
c     #############################################################
c     ##                                                         ##
c     ##  subroutine checkkernel -- check recursion kernel size  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "checkkernel" checks and updates the recursion kernel size
c     if dudl falls outside the range of maxflmda and minflmda
c
c
c      subroutine checkkernel (dudlambda)
c      use atoms
c      use mutant
c      use osrwi
c      implicit none
c      integer i,j
c      integer oldFLambdaBins,newFLambdaBins
c      integer offset
c      real*8 dudlambda
c      real*8 origDeltaG,newDeltaG,tol
c      real*8 updateFLambda
c      real*8, allocatable :: rkCopy(:,:)
c
c      tol = 1.0d-10
c      offset = 0
c
cc
cc     determine whether resize is needed
cc
c      if (dudlambda .gt. maxflmda) then
c         newFLambdaBins = FLambdaBins
c         do while (minflmda + newFLambdaBins*dFL .lt. dudlambda)
c            newFLambdaBins = newFLambdaBins + 100
c         end do
c
c      else if (dudlambda .lt. minflmda) then
c         offset = 100
c         do while (dudlambda .lt. minflmda - offset*dFL)
c            offset = offset + 100
c         end do
c         newFLambdaBins = FLambdaBins + offset
c
c      else
c         return
c      end if
c
cc
cc     save current free energy for consistency check
cc
c      origDeltaG = updateFLambda(.false.)
c      oldFLambdaBins = FLambdaBins
c
c      allocate (rkCopy(lambdaBins,oldFLambdaBins))
c      do i = 1, lambdaBins
c         do j = 1, oldFLambdaBins
c            rkCopy(i,j) = recursionKernel(i,j)
c         end do
c      end do
c
c      deallocate (recursionKernel)
c      allocate (recursionKernel(lambdaBins,newFLambdaBins))
c
cc
cc     initialize resized kernel to zero
cc
c      do i = 1, lambdaBins
c         do j = 1, newFLambdaBins
c            recursionKernel(i,j) = 0.0d0
c         end do
c      end do
c
cc
cc     copy old kernel into new kernel
cc     offset = 0 means growth on the high side
cc     offset > 0 means growth on the low side
cc
c      do i = 1, lambdaBins
c         do j = 1, oldFLambdaBins
c            recursionKernel(i,j+offset) = rkCopy(i,j)
c         end do
c      end do
c
c      FLambdaBins = newFLambdaBins
c
c      if (offset .eq. 0) then
c         maxflmda = minflmda + dFL*FLambdaBins
c         print *, 'dudlambda > maxflmda', dudlambda, maxflmda
c      else
c         minflmda = minflmda - offset*dFL
c         print *, 'dudlambda < minflmda', dudlambda, minflmda
c      end if
c
c      newDeltaG = updateFLambda(.false.)
c      if (abs(origDeltaG-newDeltaG) .gt. tol) then
c         print *, 'origDeltaG != updateFLambda'
c         print *, origDeltaG, newDeltaG
c         call fatal
c      end if
c
c      deallocate (rkCopy)
c
c      return
c      end