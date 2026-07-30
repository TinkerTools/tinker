c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine mapsublmda -- map from lambda to sublambda  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "mapsublmda" maps from the main lambda value "lmda" to the
c     electrostatic, polarization and van der Waals sublambdas
c
c
      subroutine mapsublmda (lmda)
      use dlmda
      use mutant
      implicit none
      real*8 lmda
      real*8 taper
      real*8 dtaper
      real*8 d2taper
      character*6 mode
c
c
c     map from lambda to sublambdas
c
      if (plmdamap .eq. 'EXP') then
         call sublmdaexp (lmda,plmdaexp,plambda,
     &                    dpldlmda,d2pldlmda2)
      else if (plmdamap .eq. 'INV') then
         call sublmdainvpower (lmda,plmdainvn,plmdainveps,plambda,
     &                         dpldlmda,d2pldlmda2)
      else
         mode = 'QNTPOL'
         call sublmdataper (mode,lmda,taper,dtaper,d2taper)
         plambda = 1.0d0 - taper
         dpldlmda = -dtaper
         d2pldlmda2 = -d2taper
      end if
      if (elmdamap .eq. 'EXP') then
         call sublmdaexp (lmda,elmdaexp,elambda,
     &                    deldlmda,d2eldlmda2)
      else if (elmdamap .eq. 'INV') then
         call sublmdainvpower (lmda,elmdainvn,elmdainveps,elambda,
     &                         deldlmda,d2eldlmda2)
      else
         mode = 'QNTELE'
         call sublmdataper (mode,lmda,taper,dtaper,d2taper)
         elambda = 1.0d0 - taper
         deldlmda = -dtaper
         d2eldlmda2 = -d2taper
      end if
      if (vlmdamap .eq. 'EXP') then
         call sublmdaexp (lmda,vlmdaexp,vlambda,
     &                    dvldlmda,d2vldlmda2)
      else if (vlmdamap .eq. 'INV') then
         call sublmdainvpower (lmda,vlmdainvn,vlmdainveps,vlambda,
     &                         dvldlmda,d2vldlmda2)
      else
         mode = 'QNTVDW'
         call sublmdataper (mode,lmda,taper,dtaper,d2taper)
         vlambda = 1.0d0 - taper
         dvldlmda = -dtaper
         d2vldlmda2 = -d2taper
      end if
c
c     set flags to compute polarization lambda derivative
c
      if (plmdamap .eq. 'QNT') then
         use_pol4i = (lmda .le. qntplmda1)
         use_pol4f = (lmda .ge. qntplmda0)
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine sublmdaexp -- exponential sublambda mapping  ##
c     ##                                                          ##
c     ##############################################################
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
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine sublmdainvpower -- inverse-power mapping  ##
c     ##                                                       ##
c     ###########################################################
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
c     ###############################################################
c     ##                                                           ##
c     ##  function uselmdachain -- test for an active main lambda  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "uselmdachain" is true when one of the lambda dynamics methods
c     drives a main lambda, so the sublambdas follow it through the
c     mapping and chain rule instead of being set independently
c
c
      function uselmdachain ()
      use dlmda
      implicit none
      logical uselmdachain
c
c
      uselmdachain = (use_ost .or. use_meta .or. use_ti)
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
c     #############################################################
c     ##                                                         ##
c     ##  subroutine avgstd -- average and std deviation kernel  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "avgstd" computes the average and population standard deviation
c     of the "count" samples of "list" starting at index "begin"
c
c
      subroutine avgstd (list,begin,count,avg,std)
      implicit none
      integer i
      integer begin
      integer count
      real*8 avg,std
      real*8 delta
      real*8 list(*)
c
c
c     an empty sample range has no average or deviation
c
      if (count .lt. 1) then
         avg = 0.0d0
         std = 0.0d0
         return
      end if
c
c     compute the average from the collected samples
c
      avg = 0.0d0
      do i = begin, begin+count-1
         avg = avg + list(i)
      end do
      avg = avg / dble(count)
c
c     compute the population standard deviation of the samples
c
      std = 0.0d0
      do i = begin, begin+count-1
         delta = list(i) - avg
         std = std + delta*delta
      end do
      std = sqrt(std/dble(count))
      return
      end
